from __future__ import absolute_import, print_function, unicode_literals
import sys
import os
import pandas as pd
import numpy as np
from utils import row_to_block

#import libraries for different classifiers
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn import decomposition
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_webapp')
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_analysis')
from flutype_analysis import utils
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "flutype_webapp.settings")
import django
django.setup()

class Data(object):

    def __init__(self,spots,**kwargs):
        self.spots = spots
        self.data = self._reformat_add_replica_pivot_table(**kwargs)

    @staticmethod
    def _reformat(spots):
        columns = ['raw_spot__lig_fix_batch__sid',
                  'raw_spot__lig_mob_batch__sid',
                  'intensity',
                  'raw_spot__raw_spot_collection__sid',
                  'raw_spot__row',
                  'raw_spot__column',
                  ]

        renamed_columns = ['Ligand Batch',
                           'Analyt Batch',
                           'Intensity',
                           'Collection',
                           'Row',
                           'Column']

        data = spots.values_list(*columns)
        data_numpy = np.array(data)
        spots_pd = pd.DataFrame(data_numpy,columns=renamed_columns)
        spots_pd.replace([None],[np.NaN], inplace=True)
        spots_pd["Block"] = [row_to_block(row) for row in  spots_pd["Row"]]
        return spots_pd

    @staticmethod
    def _safe_random_permuation(data, with_in):

        # on in  ["Analyt Batch","Collection","Block"]

        indexes = ["Analyt Batch", "Collection", "Block"]

        if bool(with_in):
            indexes= indexes[:indexes.index(with_in) + 1]

        safe_set = set(zip(*(data.index.get_level_values(index).values for index in indexes)))
        safe_set_indexes = [set(a) for a in safe_set]

        for i in safe_set_indexes:
            true_index =[i.issubset(set(index)) for index in data.index.values]
            data.loc[data.iloc[true_index].index] = data.iloc[true_index].apply(np.random.permutation, axis=0)
            #data.iloc[true_index].apply(np.random.permutation, axis=0)



    @staticmethod
    def _add_replica_row(spots):
        spots["Replica"] = 0
        for collection_unique in spots["Collection"].unique():
            for block_unique in spots["Block"].unique():
                for peptide_unique in spots["Ligand Batch"].unique():
                    replica = 0
                    for index in spots.index:
                        if spots["Collection"][index] == collection_unique and spots["Ligand Batch"][index] == peptide_unique and spots["Block"][index] == block_unique:
                                spots.set_value(index, "Replica", replica)
                                replica += 1


    def _reformat_add_replica_pivot_table(self,mean_on=None):

        data = self._reformat(self.spots)
        indexes = ["Analyt Batch","Collection","Block","Replica"]

        if bool(mean_on):
            indexes=indexes[:indexes.index(mean_on)+1]

        else:
            self._add_replica_row(data)

        classification_data = data.pivot_table(index=indexes,columns="Ligand Batch", values="Intensity")
        return classification_data

    def split_in_train_and_validation(self,split_on=None,ratio=None):

        #todo: add split on functionallity
        indexes = ["Analyt Batch","Collection","Block","Replica"]


        validation_data = self.data.iloc[self.data.index.get_level_values("Replica") == 0]
        train_data = self.data.drop(validation_data.index)

        return train_data , validation_data


class Analysis(object):

    def __init__(self, data):
        self.data = data
        self.classifier_names = ["Nearest Neighbors",
                 "Decision Tree", "Random Forest", "AdaBoost",
                 "Naive Bayes", "LDA"]
        self.train_data, self.validation_data = self.data.split_in_train_and_validation()

        self.y_train = pd.get_dummies(self.train_data.index.get_level_values("Analyt Batch"))
        self.y_validate = pd.get_dummies(self.validation_data.index.get_level_values("Analyt Batch"))

        self.classifiers = [
            KNeighborsClassifier(3),  # three nearest neighbors
            DecisionTreeClassifier(random_state=1),
            RandomForestClassifier(n_estimators=np.shape(self.train_data)[1], random_state=1),
            AdaBoostClassifier(),
            GaussianNB(),
            LinearDiscriminantAnalysis(n_components=5)
        ]

    @staticmethod
    def perf_measure(y_actual, y_hat):
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        y_actual = np.array(y_actual).flatten()
        y_hat = y_hat.flatten()

        for i in range(len(y_hat)):
            if y_actual[i] == y_hat[i] == 1:
                TP += 1
        for i in range(len(y_hat)):
            if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
                FP += 1
        for i in range(len(y_hat)):
            if y_actual[i] == y_hat[i] == 0:
                TN += 1
        for i in range(len(y_hat)):
            if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
                FN += 1
        return TP, FP, TN, FN


    def fit_classifiers(self):
        multi_clfs = []
        for clf in self.classifiers:
            multi_target_clf = MultiOutputClassifier(clf, n_jobs=-1)
            multi_clfs.append(multi_target_clf.fit(self.train_data,self.y_train))

        return multi_clfs

    def performance_table(self):
        perform_table = pd.DataFrame(columns=('Name', 'unweighted accuracy', 'true positive', 'false positive', 'true negative', 'false negative'))

        multi_clfs = self.fit_classifiers()

        for multi_clf, name, i in zip(multi_clfs, self.classifier_names, range(len(self.classifier_names))):
                TP, FP, TN, FN = self.perf_measure(self.y_validate, multi_clf.predict(self.validation_data))
                perform_table.loc[i] = name, 100 * multi_clf.score(self.validation_data, self.y_validate), 100 * float(TP) / (
                        TP + FP), 100 * float(FP) / (TP + FP), 100 * float(TN) / (FN + TN), 100 * float(FN) / (FN + TN)
        return perform_table


if __name__ == '__main__':
    from flutype.models import Spot
    raw_spot_collections = ["2018-01-24_E14_X31",
                            "2018-01-24_E15_X31",
                            "2018-01-24_N21_Pan",
                            "2018-01-24_N22_Cal",
                            "2018-01-24_N23_X31",
                            ]
    spots = Spot.objects.filter(raw_spot__raw_spot_collection__sid__in=raw_spot_collections)

    d=Data(spots=spots, mean_on="Block")


    ana = Analysis(d)
    print(ana.performance_table())





