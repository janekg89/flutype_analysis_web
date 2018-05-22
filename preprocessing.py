import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.interpolate import interp1d
from itertools import combinations
from spot2intensity.correlation_significance import correlated_significance
from scipy.stats import ttest_ind_from_stats


def outliers_modified_z_score(ys):
    threshold = 3.5
    median_y = np.median(ys)
    median_absolute_deviation_y = np.median([np.abs(y - median_y) for y in ys])
    modified_z_scores = [0.6745 * (y - median_y) / median_absolute_deviation_y for y in ys]
    return np.where(np.abs(modified_z_scores) > threshold)


def outlier_filtering(spots, on= "Analyte Batch"):
    frames = []
    for name, data in spots.groupby(["Ligand Batch",on]):
        out = outliers_modified_z_score(data["Intensity"])
        frames.append(data.drop(data.index[out]))
    return pd.concat(frames)


def normalize_on_ligand_batch(spots, ligand_batch = "Nenad"):
    frames = []
    for name, data in spots.groupby(["Collection"]):
        lg_int = data[data["Ligand Batch"] == ligand_batch]["Intensity"].mean()
        data["Intensity"] = data["Intensity"] / lg_int
        frames.append(data)
    return pd.concat(frames)


def normalize(spots):
    first_name = 0
    frames = []
    for name, data in spots.groupby(["Collection", "Analyte Batch"]):
        if first_name != name[1]:
            first_name = name[1]
            frames.append(data)
            c1n = data.set_index(["Ligand Batch"])
            continue
        c2n = data.set_index(["Ligand Batch"])
        x1 = np.log10(c1n["Intensity"] * c2n["Intensity"])
        y1 = np.log2(c1n["Intensity"] / c2n["Intensity"])
        lowess = sm.nonparametric.lowess(y1, x1, frac=0.33)
        f = interp1d(lowess[:, 0], lowess[:, 1], bounds_error=False)
        x_mean = x1.reset_index().groupby("Ligand Batch").agg(["mean", "std"])
        x_mean = x_mean.reset_index().rename(columns={"Intensity": "x_norm"})
        # x_norm is for normalization and is not used after normalization
        data = pd.merge(data, x_mean, how='left', left_on='Ligand Batch', right_on='Ligand Batch')
        data["Intensity"] = data["Intensity"] * (2 ** (f(data[("x_norm", "mean")])))
        frames.append(data)
    return pd.concat(frames)


def lowless_norm(spots, master_collection):
    frames = []

    c1n = spots[spots["Collection"] == master_collection]
    c1n_lf = c1n[c1n["Ligand"] == "LF"]

    for name, spots_collection in spots.groupby(["Collection"]):
        if name == master_collection:
            continue
        sp_coll_lf = spots_collection[spots_collection["Ligand"] == "LF"]
        c2n_lf = sp_coll_lf.set_index(["Ligand Batch"])
        c2n = spots_collection.set_index(["Ligand Batch"])

        x1_lf = np.log10(c1n_lf["Intensity"] * c2n_lf["Intensity"])
        x1 = np.log10(c1n["Intensity"] * c2n["Intensity"])
        y1_lf = np.log2(c1n_lf["Intensity"] / c2n_lf["Intensity"])
        lowess = sm.nonparametric.lowess(y1_lf, x1_lf, frac=0.33)
        f = interp1d(lowess[:, 0], lowess[:, 1], bounds_error=False)
        #hier comes the normlization
        x_mean = x1.reset_index().groupby("Ligand Batch").agg(["mean", "std"])
        x_mean = x_mean.reset_index().rename(columns={"Intensity": "x_norm"})
        # x_norm is for normalization and is not used after normalization
        data = pd.merge(spots_collection, x_mean, how='left', left_on='Ligand Batch', right_on='Ligand Batch')
        data["Intensity"] = data["Intensity"] * (2 ** (f(data[("x_norm", "mean")])))
        frames.append(data)
    return pd.concat(frames)

def mean_on_analyte_batch(spots):
    frames = []
    for name, data in spots.groupby(["Ligand Batch", "Analyte Batch"]):

        x = data.mean()
        x["Count"] = len(data)
        x["Intensity_std"] = data["Intensity"].std(ddof=1) / np.sqrt(len(data))
        x["Intensity_var"] = data["Intensity"].var()
        x["Intensity_rsd"] = data["Intensity"].std() / data["Intensity"].mean()
        x.name = name
        frames.append(x)
    mean_spots = pd.concat(frames, axis=1)
    return mean_spots.transpose().reset_index().rename(columns={"level_0":"Ligand Batch","level_1":"Analyte Batch"})


def mean_on_collection(spots):
    frames = []
    for name, data in spots.groupby(["Ligand Batch", "Collection","Study","Analyte Batch"]):

        x = data.mean()
        x["Count"] = len(data)
        x["Intensity_std"] = data["Intensity"].std(ddof=1) / np.sqrt(len(data))
        x["Intensity_var"] = data["Intensity"].var()
        x["Intensity_rsd"] = data["Intensity"].std() / data["Intensity"].mean()
        x.name = name
        frames.append(x)
    mean_spots = pd.concat(frames, axis=1)
    return mean_spots.transpose().reset_index().rename(columns={"level_0":"Ligand Batch","level_1":"Collection","level_2":"Study","level_3":"Analyte Batch"})


def mean_on_ligand_batch(spots):
    frames = {}
    for cn, d in spots.groupby(["Collection", "Ligand Batch"]):
        x = d.mean()
        x["Count"] = len(d)
        x["Intensity_std"] = d["Intensity"].std(ddof=1) / np.sqrt(len(d))
        x["Intensity_var"] = d["Intensity"].var()
        x["Intensity_rsd"] = d["Intensity"].std() / d["Intensity"].mean()
        frames[cn] = x
    return  pd.concat(frames, axis=1).transpose().reset_index().rename(columns={"level_1":"Ligand Batch","level_0":"Collection"})


def ligand_batch_significance(spots):
    spots_grouped = combinations(spots.groupby("Analyte Batch"), 2)
    frames = []
    for (name1, spots1), (name2, spots2) in spots_grouped:

        # get intersection of unique  ligand batches of both groups
        ligand_batches = set(spots1["Ligand Batch"].unique()).intersection(set(spots2["Ligand Batch"].unique()))

        for ligand_batch in ligand_batches:
            spots1_this_lb = spots1.loc[spots1["Ligand Batch"] == ligand_batch]
            spots2_this_lb = spots2.loc[spots2["Ligand Batch"] == ligand_batch]

            v1_i = spots1_this_lb["Intensity"].iloc[0]
            v2_i = spots2_this_lb["Intensity"].iloc[0]

            v1_err = spots1_this_lb["Intensity_std"].iloc[0]
            v2_err = spots2_this_lb["Intensity_std"].iloc[0]

            v1_var = spots1_this_lb["Intensity_var"].iloc[0]
            v2_var = spots2_this_lb["Intensity_var"].iloc[0]

            v1_count = spots1_this_lb["Count"].iloc[0]
            v2_count = spots2_this_lb["Count"].iloc[0]

            #where_max, sig_max = correlated_significance(v1_i, v2_i,num1, v1_var, v2_var,num2)
            _, sig_max = ttest_ind_from_stats(v1_i,v1_err,v1_count,v2_i,v2_err,v2_count)

            fr = pd.Series(
                [ligand_batch, (name1, name2), sig_max, v1_i, v2_i, v1_err, v2_err,v1_count,v2_count],
                index=["Ligand Batch", "Analyte Batches", "Significance", "V1_I", "V2_I", "V1_Err", "V2_Err","V1_Count", "V2_Count"])
            frames.append(fr)

    return pd.concat(frames, axis=1).transpose().dropna()

def my_aggs(x, max_spots):

    names = {
        "Intensity_rsd": np.nanmean(x["Intensity_rsd"]),
        'Count': x['Count'].sum()/max_spots[x["Analyte Batch"].unique()[0]],
        'Pixel_rsd': np.nanmean(x['Std']/x["Intensity"]), }


    return pd.Series(names)






