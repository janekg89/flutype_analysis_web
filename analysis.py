from __future__ import absolute_import, print_function, unicode_literals
import sys
import os
import pandas as pd
import numpy as np
from utils import row_to_block
import copy

import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
from collections import namedtuple

#import libraries for different classifiers
from sklearn.multioutput import MultiOutputClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

#my webapp import
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_webapp')
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_analysis')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "flutype_webapp.settings")
import django
django.setup()

class Data(object):

    def __init__(self,spots ,**kwargs):
        self.spots = spots
        #self.data = self._reformat_add_replica_pivot_table(**kwargs)
        self.kwargs = kwargs


    @staticmethod
    def _reformat(spots):
        columns = ['raw_spot__lig_fix_batch__sid',
                  'raw_spot__lig_fix_batch__concentration',
                  'raw_spot__lig_fix_batch__ligand__sid',
                  'raw_spot__lig_mob_batch__sid',
                  'raw_spot__lig_mob_batch__ligand__sid',
                  'intensity',
                  'raw_spot__raw_spot_collection__sid',
                  'raw_spot__row',
                  'raw_spot__column',
                  ]

        renamed_columns = ['Ligand Batch',
                           'Ligand Batch Concentration',
                           'Ligand',
                           'Analyte Batch',
                           'Analyte',
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
    def _safe_random_permutation(data, with_in):

        indexes = ["Analyte Batch", "Collection", "Block"]

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
        unique_instances = {}

        spots_with_replica = copy.deepcopy(spots)
        spots_with_replica["Replica"] = 0
        for i, spot in spots_with_replica[["Collection","Block","Analyte Batch","Ligand Batch"]].iterrows():

            key = tuple(spot.values)

            value = unique_instances.get(key,0)
            spots_with_replica.set_value(i,"Replica",value + 1)
            unique_instances[key] = value + 1

        return spots_with_replica


    @staticmethod
    def _append_con_diff_features(spots):
        """
        This function creates new features from existing features.
        Therefore, intensity values for the same peptides but different concentrations are subtracted.
        the return function gives the new features. If inplace True a new collumns called Constructed is added.
        For the constructed features it is set to true else false.
        :param spots:
        :return:
        """
        assert "Ligand Batch Concentration" in spots.columns, 'You have no Ligand Batch Concentration defined!'
        assert "Replica" in spots.columns, 'You have not added the Replica row!'

        spots["Constructed Feature"] = False

        for ligand_index, ligand_data in spots.groupby(["Collection","Block","Replica", "Ligand"]):

            new_feature = copy.deepcopy(ligand_data.iterrows().next()[1])

            con = ligand_data["Ligand Batch Concentration"]

            lig_low_list = ligand_data[con == con.min()]
            lig_high_list = ligand_data[con == con.max()]
            assert all([len(lig_low_list.index) == 1,len(lig_high_list.index) == 1]), "Same Concentration on with in one " \
                                                                            "Replica and Ligand exist!"
            lig_low = lig_low_list.iterrows().next()[1]
            lig_high = lig_high_list.iterrows().next()[1]

            new_feature["Intensity"] = lig_high["Intensity"] - lig_low["Intensity"]
            new_feature["Ligand Batch"] = "{}-{}".format(lig_high["Ligand Batch"],lig_low["Ligand Batch"])
            new_feature["Ligand Batch Concentration"] = np.NaN
            new_feature["Row"] = np.NaN
            new_feature["Column"] = np.NaN
            new_feature["Constructed Feature"] = True
            spots = spots.append(new_feature, ignore_index=True)

        return spots


    @staticmethod
    def _append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection"):
        spots["Normalized Feature"] = False


        constructed_data = spots[spots[which_features] == True]
        normalized_data = copy.deepcopy(constructed_data)

        for name, group in  constructed_data.groupby(with_in):
            mu, sigma = norm.fit(group["Intensity"])
            #rescalling
            datos = (group["Intensity"] - mu) / sigma
            normalized_data.set_value(group.index, "Intensity", datos)

        normalized_data["Normalized Feature"] = True
        normalized_data.reset_index(inplace=True)
        normalized_data["Ligand Batch"] = normalized_data["Ligand Batch"]+"norm"

        spots = spots.append(normalized_data, ignore_index=True)

        return spots

    @staticmethod
    def _pivot_table(spots,mean_on=None):


        indexes = ["Analyte Batch", "Collection", "Block", "Replica"]

        if bool(mean_on):
                    indexes=indexes[:indexes.index(mean_on)+1]

        for index in indexes:
            assert index in spots.columns,"Have need to create the row {}".format(index)

        data = spots.pivot_table(index=indexes,columns="Ligand Batch", values="Intensity")
        data.reset_index(inplace=True)
        return data

    def _reformat_add_replica_pivot_table(self,mean_on=None, reformat=True):

        if reformat:
            data = self._reformat(self.spots)
        else:
            data=self.spots
        indexes = ["Analyte Batch","Collection","Block","Replica"]

        if bool(mean_on):
            indexes=indexes[:indexes.index(mean_on)+1]

        else:
            self._add_replica_row(data)

        classification_data = data.pivot_table(index=indexes,columns="Ligand Batch", values="Intensity")
        classification_data.reset_index(inplace=True)
        return classification_data

    def split_in_train_and_validation(self,split_on=None,ratio=None):

        #todo: add split on functionallity
        indexes = ["Analyte Batch","Collection","Block","Replica"]


        validation_data = self.data.iloc[self.data.index.get_level_values(self.kwargs.get("mean_on", "Replica")) == 1]
        train_data = self.data.drop(validation_data.index)

        return train_data , validation_data




class Analysis(object):

    def __init__(self, data):
        self.data = data
        self.classifier_names = ["Nearest Neighbors",
                 "Decision Tree", "Random Forest", "AdaBoost",
                 "Naive Bayes", "LDA"]
        self.train_data, self.validation_data = self.data.split_in_train_and_validation()

        self.y_train = pd.get_dummies(self.train_data.index.get_level_values("Analyte Batch"))
        self.y_validate = pd.get_dummies(self.validation_data.index.get_level_values("Analyte Batch"))

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
        perform_table = pd.DataFrame(columns=('Name',
                                              'unweighted accuracy',
                                              'true positive',
                                              'false positive',
                                              'true negative',
                                              'false negative'))

        multi_clfs = self.fit_classifiers()

        for multi_clf, name, i in zip(multi_clfs, self.classifier_names, range(len(self.classifier_names))):

                TP, FP, TN, FN = self.perf_measure(self.y_validate, multi_clf.predict(self.validation_data))
                unweighted_accuracy = 100 * multi_clf.score(self.validation_data, self.y_validate)
                true_positive = 100 * np.divide(np.float64(TP) , np.float64(TP + FP))
                false_positive = 100 * np.divide(np.float64(FP) , np.float64(TP + FP))
                true_negative = 100 * np.divide(np.float64(TN) , np.float64(FN + TN))
                false_negative = 100 * np.divide(np.float64(FN) , np.float64(FN + TN))

                perform_table.loc[i] = name, \
                                       unweighted_accuracy, \
                                       true_positive, \
                                       false_positive, \
                                       true_negative, \
                                       false_negative
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





