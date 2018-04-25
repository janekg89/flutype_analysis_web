

from __future__ import absolute_import, print_function, unicode_literals
import sys
import os

#my webapp import
sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_webapp')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "flutype_webapp.settings")
import django
django.setup()

import pandas as pd
import numpy as np
import itertools
import copy
from utils import checkEqual
from preprocessing import outlier_filtering, normalize_on_ligand_batch, mean_on_analyte_batch,ligand_batch_significance
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import Imputer


class Data(object):
    def __init__(self, spots_dj = None, test=False,impute=False,spots_pd = None):
        if spots_dj is not None:
            self.spots_pd = self._reformat(spots_dj)
        else:
            assert spots_pd is not None, "Provide spots_pd or spots_dj."
            self.spots_pd = spots_pd
        self.train_test_combinations = None
        if test:
            self.train_test_combinations = self.combination_split_train_test()
        if impute:
            self.x = self._pivot_table(self.spots_pd)
        else:
            x = self._pivot_table(self.spots_pd)
            imputer = Imputer()
            imputed_x = imputer.fit_transform(x)
            self.x = pd.DataFrame(imputed_x, index=x.index, columns=x.columns)

        self.y_names = self.x.index.get_level_values("Analyte Batch").values
        self.y = self.y()
        self.collections = self.x.index.get_level_values("Collection")
        self.y_map = {i :x for i,x in enumerate(self.y.columns.tolist())}
        self.col_y_mapping = dict(zip(self.x.index.get_level_values("Collection").values,self.y_names))



    @staticmethod
    def _reformat(spots):

        columns = ['raw_spot__lig_fix_batch__sid',
                  'raw_spot__lig_fix_batch__concentration',
                  'raw_spot__lig_fix_batch__ligand__sid',
                  'raw_spot__lig_mob_batch__sid',
                  'raw_spot__lig_mob_batch__ligand__sid',
                  'intensity',
                  'std',
                  'circle_quality',
                  'raw_spot__raw_spot_collection__sid',
                  'raw_spot__raw_spot_collection__studies__sid',
                   'raw_spot__row',
                  'raw_spot__column',
                   'spot_collection__sid'


                  ]

        renamed_columns = ['Ligand Batch',
                           'Ligand Batch Concentration',
                           'Ligand',
                           'Analyte Batch',
                           'Analyte',
                           'Intensity',
                           'Std',
                           'Circle Quality',
                           'Collection',
                           'Study',
                           'Row',
                           'Column',
                           'Collection Type']

        data = spots.values_list(*columns)
        #renaming columns
        spots_pd = pd.DataFrame(np.array(data),columns=renamed_columns)
        spots_pd.replace([None], [np.NaN], inplace=True)
        spots_pd["Intensity"] = spots_pd["Intensity"].astype("float")
        return Data.add_replica(spots_pd)

    def reset_index(self):
        return Data(spots_pd=self.spots_pd.reset_index())

    def x_pca_fit(self):
        n_components = 3
        pca = PCA(n_components=n_components)
        return pca.fit(self.x)

    def x_lda_fit(self):
        n_components = 3
        clf = LinearDiscriminantAnalysis(n_components=n_components)
        return clf.fit(self.x, self.y.values.argmax(axis=1))





    @staticmethod
    def replica(spots_pd):
        spots_pd["Replica"] = range(len(spots_pd))
        return spots_pd

    @staticmethod
    def add_replica(spots_pd):
        return spots_pd.groupby(["Collection", "Ligand Batch"]).apply(Data.replica)

    def subset_collection(self, collections):
        return self._generic_subset(collections,"Collection")

    def subset_ligand(self,ligand):
        return self._generic_subset(ligand,"Ligand")

    def subset_ligand_batches(self,ligand_batches):
        return self._generic_subset(ligand_batches,"Ligand Batch")

    def _generic_subset(self,generic_list, on):
        return self.__class__(spots_pd=self.spots_pd[self.spots_pd.applymap(lambda x: x in list(generic_list))[on]])

    def collection_in_virbatch(self):
        dict_col_vir_batch = {}
        for vir_batch_name , d in self.spots_pd.groupby(["Analyte Batch"]):
            dict_col_vir_batch[vir_batch_name] = list(d["Collection"].unique())
        return dict_col_vir_batch

    def combination_split_train_test(self):
        cs = list(self.spots_pd["Collection"].unique())
        test = list(itertools.product(*self.collection_in_virbatch().values()))
        comb = pd.DataFrame(columns=["Test", "Train"])
        for i, row in enumerate(test):
            comb.set_value(i, "Test", row)
            comb.set_value(i, "Train", tuple(set(cs) - set(row)))
        return comb

    def norm(self, by="Nenad"):
        return Data(spots_pd=normalize_on_ligand_batch(self.spots_pd,ligand_batch=by))

    def outlier_filtering(self):
        return Data(spots_pd= outlier_filtering(self.spots_pd,on="Collection"))

    def mean_on_analyte_batch(self):
        return mean_on_analyte_batch(self.spots_pd)

    def ligand_batch_significance(self):
        return ligand_batch_significance(self.mean_on_analyte_batch())





    def y(self):
        return pd.get_dummies(self.y_names)


    def sample_on_analyte_batch(self,sample_size):
        return self._generic_sample("Analyte Batch",sample_size)


    def sample_on_collection(self, sample_size):
        return self._generic_sample("Collection",sample_size)

    def _generic_sample(self, on , sample_size):
        frames = []
        for name, d in self.spots_pd.groupby([on, "Ligand Batch"]):
            this_d = d.sample(sample_size, replace=True)
            this_d["Replica"] = range(sample_size)
            frames.append(this_d)
        return self.__class__(spots_pd=pd.concat(frames))


    def clean(self):
        """
        reduce to maximal common set of features
        :return:self
        """

        return self.subset_ligand_batches(self.dense_ligand_batches())


    def spare_ligand_batches(self):
        return self.x.columns[self.x.isnull().any()].tolist()

    def dense_ligand_batches(self):
        return set(self.x.columns)-set(self.spare_ligand_batches())

    def impute(self):
        imputer = Imputer()
        imputer.fit_transform(self.x)

        #return Data(spots_pd=normalize_on_ligand_batch(self.spots_pd,ligand_batch=by))

        pass



    @staticmethod
    def _pivot_table(spots_pd):
        indexes = ["Analyte Batch", "Collection", "Replica"]
        data = spots_pd.pivot_table(index=indexes,columns="Ligand Batch", values="Intensity")
        return data

    @staticmethod
    def _norm(c_data,by):
        c_data["Intensity"] = c_data["Intensity"]/c_data[c_data["Ligand Batch"] == by]["Intensity"].mean()
        return c_data


class Analysis(object):

    def __init__(self, data, **kwargs ):

        self.train_test = kwargs.get('train_test', data.train_test_combinations)
        self.data = data

        self.classifier_names = [#"Nearest Neighbors",
                                 #"Decision Tree",
                                 # "Random Forest",
                                 # "AdaBoost",
                                 #"Gaussian NB",
                                 #"LDA",
                                "LogisticRegression",
                                 ]

        self.classifiers = [
            #KNeighborsClassifier(3),  # three nearest neighbors
            #DecisionTreeClassifier(random_state=1),
            #RandomForestClassifier(n_estimators=np.shape(self.train_data)[1], random_state=1),
            #AdaBoostClassifier(),
            #GaussianNB(),
            #GaussianNB(),
            #LinearDiscriminantAnalysis(n_components=5),
            LogisticRegression(multi_class ="multinomial",solver='lbfgs'),
        ]


    def calculate_all(self):
         frames = []

         for c_name, c in zip(self.classifier_names,self.classifiers):
             print("*"*5+c_name+"*"*5)
             d = copy.deepcopy(self.train_test)

             d["Model"] = c
             d["Classifier Name"] = c_name
             print("*"*5+"Fit Models"+"*"*5)
             d["Model"] = d.apply(self._fit, axis=1, args=(self.data,))
             print("*"*5+"Predict"+"*"*5)
             d["Predictions"] = d.apply(self._predict , axis=1, args =(self.data,))
             #print("*"*5+"Confusion Matrix"+"*"*5)

             #d["Confusion Matrix"] = d.apply(self._confusion_matrix, axis=1, args=(self.data,))
             #print("*"*5+"Score"+"*"*5)

             #d["Score"] = d.apply(self._score, axis=1, args=(self.data,))
             #print("*"*5+"Score by Collection"+"*"*5)

             #d["Score by Collection"] = d.apply(self._score_by_collection,axis=1, args=(self.data,))
             #print("*"*5+"Majority Score by Collection"+"*"*5)

             #d["Majority Score by Collection"] = d.apply(self._majority_score_by_collection, axis=1, args=(self.data,))
             #print("*"*5+"True False"+"*"*5)

             d["True False"] =d.apply(self._predicted_count_by_analyte_batch, axis=1, args=(self.data,))

             #print("*"*5+"Confusion Matrix by Collection"+"*"*5)

             #d["Confusion Matrix by Collection"] = d.apply(self._confusion_matrix_by_collection, axis=1, args=(self.data,))

             frames.append(d)
         self.result = pd.concat(frames, ignore_index=True)

    def complete_information(self):
        return pd.concat([pd.DataFrame(d) for d in self.result["True False"]])



    def subset_collection(self,collections):
        return self.data.subset_collection(collections)

    def box_plot(self, **kwargs):
        assert hasattr(self, "result"), "First run .calcualte_all()"
        self.result.boxplot(column="Score", by="Classifier Name")

    def plot_confusion_matrixes(self,normalize=False, **kwargs):
        pass


    def confusion_matrix(self):
        assert hasattr(self, "result"), "First run .calcualte_all()"
        frames = []
        for n , d in self.result.groupby("Classifier Name"):
             d =  self.result["Confusion Matrix"].apply(lambda x: pd.DataFrame(list(x), index=self.data.y.columns, columns=self.data.y.columns)).sum().apply(self.norm,axis=1)
             d["Classifier Name"] = n
             frames.append(d)

        return pd.concat(frames,levels="Classifier Name").reset_index().rename(columns={"index":"Analyte Batch"}).set_index(["Classifier Name","Analyte Batch"])


    def score_by_collection(self,stats=True):
        score = self.result["Score by Collection"].apply(lambda x: pd.Series(x))
        score.index = self.result["Classifier Name"]
        frames = {}
        for n, d in score.groupby(score.index):
            d_classifier = d.apply(self.stats)
            d_classifier.name = n
            frames[n] = d_classifier
        if stats:
            return pd.concat(frames)

        return score


    def score(self):
        df = self.score_by_collection(False).set_index(self.result["Classifier Name"]).transpose()
        df["Analyte Batch"] = [self.data.col_y_mapping[x] for x in df.index.values]
        df = pd.DataFrame(df.set_index(["Analyte Batch", df.index]).stack(), columns=["Score"])
        df.index.names = ['Analyte Batch', "Collection", "Classifier Name"]
        return df

    def score_mean(self):
        return self.score().groupby(["Classifier Name","Analyte Batch"]).mean()

    def score_std(self):
        return self.score().groupby(["Classifier Name","Analyte Batch"]).std().rename(columns={"Score":"Score Std"})

    def majority_score(self):
        df = self.result["Majority Score by Collection"].apply(pd.Series).set_index(self.result["Classifier Name"]).transpose()
        df["Analyte Batch"] = [self.data.col_y_mapping[x] for x in df.index.values]
        df = pd.DataFrame(df.set_index(["Analyte Batch",df.index]).stack() , columns=["Majority Score"])
        df.index.names = ['Analyte Batch', "Collection", "Classifier Name"]
        return df

    def majority_score_mean(self):
         return self.majority_score().groupby(["Classifier Name", "Analyte Batch"]).mean()


    def all_score(self):
        return pd.concat([self.score_mean(),self.score_std(),self.majority_score_mean()], axis=1)


    @staticmethod
    def _fit(train_test_row, data):
        train_data = data.subset_collection(train_test_row["Train"])
        c = train_test_row["Model"]
        c.fit(train_data.x, train_data.y_names)
        return c

    @staticmethod
    def _predict(train_test_row, data):
        test_data = data.subset_collection(train_test_row["Test"])
        c = train_test_row["Model"]
        return tuple(c.predict(test_data.x))


    @staticmethod
    def _confusion_matrix(train_test_row, data):
        test_data = data.subset_collection(train_test_row["Test"])
        return tuple(confusion_matrix(test_data.y_names,train_test_row["Predictions"],labels=test_data.y.columns))

    @staticmethod
    def _score(train_test_row, data):
        test_data = data.subset_collection(train_test_row["Test"])
        return train_test_row["Model"].score(test_data.x,test_data.y_names)

    @staticmethod
    def _score_by_collection(train_test_row, data):
        frames = {}
        test_data = data.subset_collection(train_test_row["Test"])
        for c_name  in test_data.spots_pd["Collection"].unique():
            c_data = test_data.subset_collection([c_name])
            frames[c_name] = train_test_row["Model"].score( c_data.x,  c_data.y_names)
        return frames

    @staticmethod
    def _confusion_matrix_by_collection(train_test_row, data):
        test_data = data.subset_collection(train_test_row["Test"])
        frames = {}
        test_compersion = pd.DataFrame(test_data.y_names, index= test_data.x.index, columns=["y"])
        test_compersion["Predictions"] = train_test_row["Predictions"]
        for c_name , c_data in test_compersion.groupby(test_compersion.index.get_level_values("Collection")):


            frames[c_name] = confusion_matrix(c_data.y.values,c_data["Predictions"].values, labels=test_data.y.columns)

        return frames

    @staticmethod
    def _majority_score_by_collection(train_test_row,data):
        test_data = data.subset_collection(train_test_row["Test"])
        frames = {}
        test_compersion = pd.DataFrame(test_data.y_names, index=test_data.x.index, columns=["y"])
        test_compersion["Predictions"] = train_test_row["Predictions"]
        for c_name, c_data in test_compersion.groupby(test_compersion.index.get_level_values("Collection")):
             a_batches = c_data.index.get_level_values("Analyte Batch").values
             assert checkEqual(a_batches), "Collection has multiple Analyte Batches!"
             a_batch = a_batches[0]
             count_analyte_batch = pd.Series(c_data["Predictions"]).value_counts()
             if len(np.argwhere(count_analyte_batch == np.amax(count_analyte_batch)))>1:
                 frames[c_name] = 0
             else:
                 frames[c_name] = int(count_analyte_batch.idxmax() == a_batch)

        return frames

    @staticmethod
    def _predicted_count_by_analyte_batch(train_test_row, data):
        test_data = data.subset_collection(train_test_row["Test"])
        test_comperison = pd.DataFrame(test_data.y_names, index=test_data.x.index, columns=["y"])
        test_comperison["Predictions"] = train_test_row["Predictions"]
        test_comperison["Test"] = [(train_test_row["Test"])]*len(test_comperison["Predictions"])

        test_comperison["TrueFalse"] = (test_comperison.index.get_level_values("Analyte Batch") == test_comperison["Predictions"]).astype(int)

        return test_comperison.to_dict()


    @staticmethod
    def norm(data):
        return data / data.sum()

    @staticmethod
    def stats(column):
        c = pd.Series()
        c["mean"] = column.mean()
        c["std"] = column.std()
        c["count"] = column.count()
        return c



















