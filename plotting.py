from __future__ import absolute_import, print_function, unicode_literals
import sys
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import analysis as a
import itertools
import utils
import copy


sys.path.append('/home/janekg89/Develop/Pycharm_Projects/flutype_webapp')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "flutype_webapp.settings")
import django
django.setup()
from flutype.models import Spot

import matplotlib.lines as mlines
from collections import OrderedDict
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def plot_block_correlation_mean_linked_concentrations(data,blocks,name, max_unzoomed):

    data1 = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand']).mean()
    data1["Std"] = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand'])['Intensity'].std() / np.sqrt(3.0)
    data=data1
    data.dropna(inplace=True)
    data.reset_index(inplace=True)

    markers = ["x","o","v","<",">","^"]

    concentrations = data["Ligand Batch Concentration"].unique()
    marker = OrderedDict(zip(markers,concentrations))


    handels_marker = [mlines.Line2D([], [], color='k', marker=m, linestyle='None', markersize=10, label=marker[m]) for m
                      in marker]
    n = len(data.groupby('Ligand').groups)
    color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))

    handels = []

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    axes = (ax1, ax2)

    for ligand_name, ligand_data in data.groupby('Ligand'):
        c = next(color)
        patch = mpatches.Patch(color=c, label=ligand_name)
        handels.append(patch)

        for ax in axes:
            marker_it = iter(marker)
            x = ligand_data["Intensity"].loc[data['Block'] == blocks[0]]
            y = ligand_data["Intensity"].loc[data['Block'] == blocks[1]]
            ax.plot(x, y,marker=None, alpha=.85, c=c, linewidth=2)

            for lb_name, lb_data in ligand_data.groupby('Ligand Batch Concentration'):
                m = next(marker_it)

                x = lb_data["Intensity"].loc[data['Block'] == blocks[0]]
                y = lb_data["Intensity"].loc[data['Block'] == blocks[1]]
                x_std = lb_data["Std"].loc[data['Block'] == blocks[0]]
                y_std = lb_data["Std"].loc[data['Block'] == blocks[1]]
                ax.errorbar(x, y, xerr=x_std, yerr=y_std, alpha=.4, capthick=2, c="black", linewidth=0.5)
                ax.scatter(x, y, alpha=.9, marker=m, c=c)

    # unzoomed axes
    max_i = max_unzoomed
    ax1.plot([0, max_i], [0, max_i], 'k--')

    # zoomed axes
    max_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].max()
    min_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].min()
    ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

    for ax in (ax1,):
        l1 = ax.legend(handles=handels_marker, bbox_to_anchor=(1.0, 0.13))
        l2 = ax.legend(handles=handels)
        ax.add_artist(l1)
        ax.add_artist(l2)

    for ax in axes:
        ax.set_xlabel('Block {}'.format(blocks[0]))
        ax.set_ylabel('Block {}'.format(blocks[1]))
        ax.set_title('Blocks {} of {}'.format(blocks, name))

    return fig


def plot_block_correlation(data,blocks,name,max_unzoomed):

    data.dropna(inplace=True)

    markers = ["x","o","v","<",">","^"]
    concentrations = data["Ligand Batch Concentration"].unique()
    marker = OrderedDict(zip(markers,concentrations))
    handels_marker = [mlines.Line2D([], [], color='k', marker=m, linestyle='None', markersize=10, label=marker[m]) for m
                      in marker]
    n = len(data.groupby('Ligand').groups)
    color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
    handels = []

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    axes = (ax1, ax2)

    for ligand_name, ligand_data in data.groupby('Ligand'):
        c = next(color)
        patch = mpatches.Patch(color=c, label=ligand_name)
        handels.append(patch)

        for ax in axes:
            marker_it = iter(marker)
            for lb_name, lb_data in ligand_data.groupby('Ligand Batch'):
                m = next(marker_it)

                x = lb_data["Intensity"].loc[data['Block'] == blocks[0]]
                y = lb_data["Intensity"].loc[data['Block'] == blocks[1]]
                ax.scatter(x, y, alpha=.9, marker=m, c=c)

    # unzoomed axes
    max_i = max_unzoomed
    ax1.plot([0, max_i], [0, max_i], 'k--')

    # zoomed axes
    max_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].max()
    min_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].min()
    ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

    for ax in (ax1,):
        l1 = ax.legend(handles=handels_marker, bbox_to_anchor=(1.0, 0.13))
        l2 = ax.legend(handles=handels)
        ax.add_artist(l1)
        ax.add_artist(l2)

    for ax in axes:
        ax.set_xlabel('Block {}'.format(blocks[0]))
        ax.set_ylabel('Block {}'.format(blocks[1]))
        ax.set_title('Blocks {} of {}'.format(blocks, name))

    return fig

def plot_array_block_correlation(data,col_pair,study):

    data1 = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand']).mean()
    data1["Std"] = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand'])['Intensity'].std() / np.sqrt(3.0)
    data = data1
    data.reset_index(inplace=True)
    marker = OrderedDict([('x', 0.25), ('o', 0.5), ('v', 1.0)])

    handels_marker = [mlines.Line2D([], [], color='k', marker=m, linestyle='None', markersize=10, label=marker[m]) for m
                      in marker]

    for block_name, block_data in data.groupby('Block'):
        n = len(block_data.groupby('Ligand').groups)

        color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
        handels = []

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        axes = (ax1, ax2)

        for ligand_name, ligand_data in block_data.groupby('Ligand'):
            c = next(color)
            patch = mpatches.Patch(color=c, label=ligand_name)
            handels.append(patch)

            for ax in axes:
                marker_it = iter(marker)
                x = ligand_data[ligand_data["Collection"] == col_pair[0]]["Intensity"]
                y = ligand_data[ligand_data["Collection"] == col_pair[1]]["Intensity"]
                ax.plot(x, y, marker=None, alpha=.85, c=c, linewidth=2)

                for lb_name, lb_data in ligand_data.groupby('Ligand Batch Concentration'):
                    m = next(marker_it)

                    x = lb_data[lb_data["Collection"] == col_pair[0]]["Intensity"]
                    y = lb_data[lb_data["Collection"] == col_pair[1]]["Intensity"]
                    x_std = lb_data[lb_data["Collection"] == col_pair[0]]["Std"]
                    y_std = lb_data[lb_data["Collection"] == col_pair[1]]["Std"]
                    ax.errorbar(x, y, xerr=x_std, yerr=y_std, alpha=.4, capthick=2, c="black", linewidth=0.5)
                    ax.scatter(x, y, alpha=.9, marker=m, c=c)

        # unzoomed axes
        max_i = data["Intensity"].max()
        min_i = data["Intensity"].min()
        ax1.plot([0, max_i], [0, max_i], 'k--')

        # zoomed axes
        max_i = block_data[[collection in col_pair for collection in block_data.Collection]]["Intensity"].max()
        min_i = block_data[[collection in col_pair for collection in block_data.Collection]]["Intensity"].min()
        ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

        for ax in (ax1,):
            l1 = ax.legend(handles=handels_marker, bbox_to_anchor=(1.0, 0.13))
            l2 = ax.legend(handles=handels)
            ax.add_artist(l1)
            ax.add_artist(l2)

        for ax in axes:
            ax.set_xlabel('{}'.format(col_pair[0]))
            ax.set_ylabel('{}'.format(col_pair[1]))
            ax.set_title('Block {} of {}'.format(block_name, col_pair))

        file_path = 'results/{}/block_correlations_new/{}_block_{}'.format(study, col_pair, block_name)
        utils.ensure_dir(file_path)
        print(file_path)

        plt.savefig(file_path)

def plot_array_correlation(data,col_pair,study,append_name=""):

    data1 = data.groupby(['Collection', 'Ligand Batch', 'Ligand']).mean()
    data1["Std"] = data.groupby(['Collection', 'Ligand Batch', 'Ligand'])['Intensity'].std() / np.sqrt(3.0)
    data = data1
    data.reset_index(inplace=True)
    marker = OrderedDict([('x', 0.25), ('o', 0.5), ('v', 1.0)])

    handels_marker = [mlines.Line2D([], [], color='k', marker=m, linestyle='None', markersize=10, label=marker[m]) for m
                      in marker]


    n = len(data.groupby('Ligand').groups)

    color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
    handels = []

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    axes = (ax1, ax2)

    for ligand_name, ligand_data in data.groupby('Ligand'):
        c = next(color)
        patch = mpatches.Patch(color=c, label=ligand_name)
        handels.append(patch)

        for ax in axes:
            marker_it = iter(marker)
            x = ligand_data[ligand_data["Collection"] == col_pair[0]]["Intensity"]
            y = ligand_data[ligand_data["Collection"] == col_pair[1]]["Intensity"]
            ax.plot(x, y, marker=None, alpha=.85, c=c, linewidth=2)

            for lb_name, lb_data in ligand_data.groupby('Ligand Batch Concentration'):
                m = next(marker_it)

                x = lb_data[lb_data["Collection"] == col_pair[0]]["Intensity"]
                y = lb_data[lb_data["Collection"] == col_pair[1]]["Intensity"]
                x_std = lb_data[lb_data["Collection"] == col_pair[0]]["Std"]
                y_std = lb_data[lb_data["Collection"] == col_pair[1]]["Std"]
                ax.errorbar(x, y, xerr=x_std, yerr=y_std, alpha=.4, capthick=2, c="black", linewidth=0.5)
                ax.scatter(x, y, alpha=.9, marker=m, c=c)

    # unzoomed axes
    max_i = data["Intensity"].max()
    min_i = data["Intensity"].min()
    ax1.plot([0, max_i], [0, max_i], 'k--')

    # zoomed axes
    max_i = data[[collection in col_pair for collection in data.Collection]]["Intensity"].max()
    min_i = data[[collection in col_pair for collection in data.Collection]]["Intensity"].min()
    ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

    for ax in (ax1,):
        l1 = ax.legend(handles=handels_marker, bbox_to_anchor=(1.0, 0.13))
        l2 = ax.legend(handles=handels)
        ax.add_artist(l1)
        ax.add_artist(l2)

    for ax in axes:
        ax.set_xlabel('{}'.format(col_pair[0]))
        ax.set_ylabel('{}'.format(col_pair[1]))
        ax.set_title('Correltation of {}'.format(col_pair))

    file_path = 'results/{}/block_correlations_new/{}_{}'.format(study, col_pair,append_name)
    print(file_path)
    utils.ensure_dir(file_path)
    plt.savefig(file_path)

def plot_array_correlation_constructed(data, col_pair, study):

        data1 = data.groupby(['Collection', 'Ligand Batch','Block', 'Ligand']).mean()
        data1["Std"] = data.groupby(['Collection', 'Ligand Batch', 'Block','Ligand'])['Intensity'].std() / np.sqrt(3.0)
        data = data1
        data.reset_index(inplace=True)
        marker = OrderedDict([('x', "Block 1"), ('o', "Block 2"), ('v','Block 3')])

        handels_marker = [mlines.Line2D([], [], color='k', marker=m, linestyle='None', markersize=10, label=marker[m])
                          for m
                          in marker]

        n = len(data.groupby('Ligand').groups)

        color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
        handels = []

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        axes = (ax1, ax2)

        for ligand_name, ligand_data in data.groupby('Ligand'):
            c = next(color)
            patch = mpatches.Patch(color=c, label=ligand_name)
            handels.append(patch)

            for ax in axes:
                marker_it = iter(marker)
                x = ligand_data[ligand_data["Collection"] == col_pair[0]]["Intensity"]
                y = ligand_data[ligand_data["Collection"] == col_pair[1]]["Intensity"]
                ax.plot(x, y, marker=None, alpha=.85, c=c, linewidth=2)

                for lb_name, lb_data in ligand_data.groupby('Ligand Batch'):

                    lb_data.sort_values(by=['Block'],inplace=True)
                    x = lb_data[lb_data["Collection"] == col_pair[0]]["Intensity"]
                    y = lb_data[lb_data["Collection"] == col_pair[1]]["Intensity"]
                    x_std = lb_data[lb_data["Collection"] == col_pair[0]]["Std"]
                    y_std = lb_data[lb_data["Collection"] == col_pair[1]]["Std"]
                    ax.errorbar(x, y, xerr=x_std, yerr=y_std, alpha=.4, capthick=2, c="black", linewidth=0.5)
                    ax.scatter(x, y, alpha=.9, c=c)

        # unzoomed axes
        max_i = data["Intensity"].max()
        min_i = data["Intensity"].min()
        ax1.plot([min_i, max_i], [min_i, max_i], 'k--')

        # zoomed axes
        max_i = data[[collection in col_pair for collection in data.Collection]]["Intensity"].max()
        min_i = data[[collection in col_pair for collection in data.Collection]]["Intensity"].min()
        ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

        for ax in (ax1,):
            l2 = ax.legend(handles=handels)
            ax.add_artist(l2)

        for ax in axes:
            ax.set_xlabel('{}'.format(col_pair[0]))
            ax.set_ylabel('{}'.format(col_pair[1]))
            ax.set_title('Constructed {}'.format(col_pair))

        file_path = 'results/{}/block_correlations_new/{}_constructed_normalized'.format(study, col_pair)
        utils.ensure_dir(file_path)
        print(file_path)
        plt.savefig(file_path)

def plot_array_block_correlation_constructed(data,col_pair,study):

    data1 = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand']).mean()
    data1["Std"] = data.groupby(['Collection', 'Block', 'Ligand Batch', 'Ligand'])['Intensity'].std() / np.sqrt(3.0)
    data = data1
    data.reset_index(inplace=True)

    for block_name, block_data in data.groupby('Block'):
        n = len(block_data.groupby('Ligand').groups)

        color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
        handels = []

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
        axes = (ax1, ax2)

        for ligand_name, ligand_data in block_data.groupby('Ligand'):
            c = next(color)
            patch = mpatches.Patch(color=c, label=ligand_name)
            handels.append(patch)

            for ax in axes:
                x = ligand_data[ligand_data["Collection"] == col_pair[0]]["Intensity"]
                y = ligand_data[ligand_data["Collection"] == col_pair[1]]["Intensity"]
                x_std = ligand_data[ligand_data["Collection"] == col_pair[0]]["Std"]
                y_std = ligand_data[ligand_data["Collection"] == col_pair[1]]["Std"]
                ax.errorbar(x, y, xerr=x_std, yerr=y_std, alpha=.4, capthick=2, c="black", linewidth=0.5)
                ax.scatter(x,y,c=c)

        # unzoomed axes
        max_i = data["Intensity"].max()
        min_i = data["Intensity"].min()
        ax1.plot([min_i, max_i], [min_i, max_i], 'k--')

        # zoomed axes
        max_i = block_data[[collection in col_pair for collection in block_data.Collection]]["Intensity"].max()
        min_i = block_data[[collection in col_pair for collection in block_data.Collection]]["Intensity"].min()
        ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

        for ax in (ax1,):
            l2 = ax.legend(handles=handels)
            ax.add_artist(l2)

        for ax in axes:
            ax.set_xlabel('{}'.format(col_pair[0]))
            ax.set_ylabel('{}'.format(col_pair[1]))
            ax.set_title('Block {} of {}'.format(block_name, col_pair))

        file_path = 'results/{}/block_correlations_new/{}_block_{}_const_norm'.format(study, col_pair, block_name)
        print(file_path)
        utils.ensure_dir(file_path)
        plt.savefig(file_path)

def plot_block_correlation_constructed(data,blocks,name,max_unzoomed):

    data.reset_index(inplace=True)

    n = len(data.groupby('Ligand').groups)
    color = iter(plt.cm.tab20c(np.linspace(0, 1, n)))
    handels = []

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    axes = (ax1, ax2)

    for ligand_name, ligand_data in data.groupby('Ligand'):
        c = next(color)
        patch = mpatches.Patch(color=c, label=ligand_name)
        handels.append(patch)

        for ax in axes:


            x = ligand_data["Intensity"].loc[ ligand_data['Block'] == blocks[0]]
            y = ligand_data["Intensity"].loc[ ligand_data['Block'] == blocks[1]]
            ax.scatter(x, y, alpha=.9, c=c)

    # unzoomed axes
    max_i = max_unzoomed
    ax1.plot([-max_i, max_i], [-max_i, max_i], 'k--')

    # zoomed axes
    max_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].max()
    min_i = data["Intensity"].iloc[[b in blocks for b in data["Block"]]].min()
    ax2.plot([min_i, max_i], [min_i, max_i], 'k--')

    for ax in (ax1,):
        l2 = ax.legend(handles=handels)
        ax.add_artist(l2)

    for ax in axes:
        ax.set_xlabel('Block {}'.format(blocks[0]))
        ax.set_ylabel('Block {}'.format(blocks[1]))
        ax.set_title('Blocks {} of {}'.format(blocks, name))

    return fig


if __name__ == '__main__':
    from flutype.models import Spot
    raw_spot_collections = ["2018-01-24_E14_X31",
                            "2018-01-24_E15_X31",
                            "2018-01-24_N21_Pan",
                            "2018-01-24_N22_Cal",
                            "2018-01-24_N23_X31",
                            ]

    study = "2018-01-24_microarray"
    spots = Spot.objects.filter(raw_spot__raw_spot_collection__sid__in=raw_spot_collections)
    spots = spots.filter(spot_collection__sid="quant1")

    d=a.Data(spots=spots)
    spots = d._reformat(spots)
    spots = d._add_replica_row(spots)
    # spots = d._append_con_diff_features(spots)
    # spots = d._append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection")
    max_unzoomed = spots["Intensity"].max()
    #data = d._pivot_table(spots)
    blocks_sets = [[1, 2], [2, 3], [1, 3]]

    """
    for blocks in blocks_sets:
        for collection_name, collection_data in spots.groupby('Collection'):

            fig = plot_block_correlation_mean_linked_concentrations(collection_data, blocks, collection_name, max_unzoomed)
            file_path = 'results/{}/block_correlations_new/{}_blocks_{}_mean_on_block'.format(study,collection_name,blocks)
            utils.ensure_dir(file_path)
            print(file_path)
            plt.savefig(file_path)


      


    
    for blocks in blocks_sets:
        for collection_name, collection_data in spots.groupby('Collection'):
            fig = plot_block_correlation(collection_data, blocks, collection_name, max_unzoomed)
            file_path = 'results/{}/block_correlations_new/{}_blocks_{}'.format(study,collection_name,blocks)
            utils.ensure_dir(file_path)
            print(file_path)
            plt.savefig(file_path)

    
            
    collections = spots.groupby("Collection").groups
    col_pairs = itertools.combinations(collections, 2)
    for col_pair in col_pairs:
        fig = plot_array_block_correlation(spots,col_pair=col_pair,study=study)


    
    collections = spots.groupby("Collection").groups
    col_pairs = itertools.combinations(collections, 2)
    for col_pair in col_pairs:
        fig = plot_array_correlation(spots,col_pair=col_pair,study=study)
        
     

    spots = d._append_con_diff_features(spots)
    spots = d._append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection")
    spots = spots[spots["Normalized Feature"]]
    collections = spots.groupby("Collection").groups
    col_pairs = itertools.combinations(collections, 2)
    for col_pair in col_pairs:
        plot_array_correlation_constructed(spots, col_pair=col_pair, study=study)
    collections = spots.groupby("Collection").groups
    col_pairs = itertools.combinations(collections, 2)

  
    for col_pair in col_pairs:
        plot_array_block_correlation_constructed(spots, col_pair=col_pair, study=study)

    """
    

    spots = d._append_con_diff_features(spots)
    spots = d._append_normalized_features(spots, which_features="Constructed Feature", with_in="Collection")
    spots = spots[spots["Normalized Feature"]]
    max_unzoomed= spots["Intensity"].max()
    for blocks in blocks_sets:
        for collection_name, collection_data in spots.groupby('Collection'):
            fig = plot_block_correlation_constructed(collection_data, blocks, collection_name, max_unzoomed)
            file_path = 'results/{}/block_correlations_new/{}_blocks_{}_const_norm'.format(study, collection_name, blocks)
            utils.ensure_dir(file_path)
            print(file_path)
            plt.savefig(file_path)
    
    




