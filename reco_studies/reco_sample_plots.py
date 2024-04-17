#!/usr/bin/env python3

# LOAD Packages 
import uproot
import pandas as pd
import numpy as np
import awkward as ak
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#matplotlib.rcParams.update(matplotlib.rcParamsDefault)
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from matplotlib import cm, colors
import matplotlib.patches as mpatches
import h5py
import argparse
from matplotlib.backends.backend_pdf import PdfPages

from plot_utils import rasterize_plots, vectorize_plots
rasterize_plots()

def main(ttree_file, protons_only = False):

    tracks = uproot.open(ttree_file)
    df = tracks["RecoBenchmarkTree"].arrays(library="pd")

    # Bin tracks by length + add binning to dataframe
    num_length_bins = 15
    max_length = 150
    length_bins = np.linspace(0.,max_length, num_length_bins)
    reco_length_bin_nums = np.digitize(df['reco_length'], length_bins)
    df['reco_length_bin'] = reco_length_bin_nums

    # Filter by overlap between reco and true tracks
    ovlp_cut = 0.5
    df_filtered = df[df['overlap']>=ovlp_cut]

    # Set up sample labels
    if protons_only:
        df_filtered = df_filtered[df_filtered['true_pdg']==2212]
        sample = '_protons_'
        sample_label = 'Protons'
    else:
        sample = '_charged_tracks_'
        sample_label = 'Charged Tracks'

    # Group df_filtered by reco_length_bin for relevant plots
    grouped_length = df_filtered.groupby('reco_length_bin')

    print('\n----------------- File content -----------------')
    print('File:',ttree_file)
    print('Sample:', sample_label)
    print('Number of tracks:', len(df_filtered))
    print('------------------------------------------------\n')

    # Set up output file
    output_pdf_name = ttree_file.split('.root')[0]+sample+str(ovlp_cut)+'_ovlp_cut_validations.pdf'
    # put output in this directory vs. where TTree file is stored
    output_pdf_name = output_pdf_name.split('/')[-1] # !!

    # Make plots in output file
    with PdfPages(output_pdf_name) as output:

        # Plot Angles by Length for various numbers of bins shown
        # Also Plot differences in angles between true and reco tracks
        angle_names = ['angle', 'angle_x', 'angle_y', 'angle_z']
        cut_lines = [-0.5, 0.75, -0.5,-0.5]
        angle_labels = ['Track Angle w.r.t. Beam', 'Track Angle w.r.t. X-Axis', \
                                 'Track Angle w.r.t. Y-Axis', 'Track Angle w.r.t. Z-Axis']
        bins_to_show = [5, 15]
        for j in range(len(angle_names)):

            # Difference Plot
            fig = plt.figure(figsize=[8,6])
            axs = fig.add_subplot(111)
            diff_counts, diff_bins, _ = axs.hist((abs(np.cos(df_filtered['true_'+angle_names[j]]))\
                                                  -abs(np.cos(df_filtered['reco_'+angle_names[j]]))), \
                                                  alpha=1.0, bins=100, range=(-1, 1), histtype='stepfilled', \
                                                  color='limegreen')
            axs.set_xlabel('Difference in True vs. ML Reco Absolute Value of Cosine of '+angle_labels[j], fontsize=12)
            axs.set_ylabel('Count [Tracks / 0.02]', fontsize=12)
            axs.set_xlim(-1,1)
            axs.set_yscale('log')
            plt.axvline(x=cut_lines[j], color='navy', linestyle='--') 
            fig.suptitle('(True - ML Reco) Difference in Absolute Value of Cosine of \n'+angle_labels[j]+' for Reconstructed '+sample_label+' Sample', fontsize=16)
            plt.tight_layout()
            output.savefig()
            plt.close()

            # Angles by Track Length plots
            for b in bins_to_show:

                plot_angle_by_length(sample_name=sample_label,angle_field=angle_names[j], \
                                     angle_label=angle_labels[j], \
                                     df_length_grouped = grouped_length, \
                                     bins_shown = b, max_bins = num_length_bins, \
                                     max_track_length = max_length)
                output.savefig()
                plt.close()



        # Plot Track Lengths
        fig = plt.figure(figsize=[8,6])
        axs = fig.add_subplot(111)
        reco_counts, reco_bins, _ = axs.hist(df_filtered['reco_length'], alpha=0.5, bins=50, \
                                            range=(0, 150), histtype='stepfilled', \
                                            label='ML Reco Track Length')
        true_counts, true_bins, _ = axs.hist(df_filtered['true_length'], alpha=0.5, bins=50, \
                                            range=(0, 150),histtype='stepfilled', \
                                            label='True Track Length')
        axs.set_xlabel('Track Length [cm]', fontsize=12)
        axs.set_ylabel('Count ['+sample_label+' / 3 cm]', fontsize=12)
        axs.set_xlim(0, 150)
        axs.legend(loc='upper right', fontsize=12)
        axs.set_title('True vs. ML Reco Track Length for Reconstructed'+sample_label+' Sample', \
                      fontsize=15)
        plt.tight_layout()
        output.savefig()
        plt.close()

        # Plot Track Lengths Difference and 2D hist
        fig = plt.figure(figsize=[8,6])
        axs = fig.add_subplot(111)
        diff_counts, diff_bins, _ = axs.hist(df_filtered['true_length']-df_filtered['reco_length'], \
                                             alpha=1.0, bins=150, range=(-150, 150),\
                                             histtype='stepfilled', color='limegreen')
        axs.set_xlabel('Difference in True vs. Reco Track Length [cm]', fontsize=12)
        axs.set_ylabel('Count ['+sample_label+' / 2 cm]', fontsize=12)
        axs.set_xlim(-130, 75)
        axs.set_yscale('log')
        plt.axvline(x=10, color='navy', linestyle='--')  # Draw a red dashed line at x=10
        plt.axvline(x=-10, color='navy', linestyle='--')  # Draw a red dashed line at x=-10
        fig.suptitle('(True - ML Reco) Track Length Difference \nfor Reconstructed '+sample_label+' Sample', fontsize=16)
        plt.tight_layout()  
        output.savefig()
        plt.close()

        vectorize_plots()
        fig, ax = plt.subplots(figsize=[8,6])
        counts, xedges, yedges, im = ax.hist2d(df_filtered['reco_length'], \
                                               df_filtered['true_length'], \
                                               bins=50, range=[[0, 150], [0, 150]], cmap='pink')
        ax.set_xlabel('ML Reco Track Length [cm]', fontsize=12)
        ax.set_ylabel('True Track Length [cm]', fontsize=12)
        ax.set_title('ML Reco vs True Match Length for Reconstructed '+sample_label+' Sample', fontsize=16)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Counts ['+sample_label+' / 3 cm]', fontsize=12)
        plt.tight_layout()
        output.savefig()
        plt.close()
        rasterize_plots()

        # Plot Track Start and End Points in XYZ
        positions_names = ['_start_', '_end_']
        positions_labels = ['Start', 'End']
        for j in range(len(positions_names)):
            pos = positions_names[j]

            fig, axs = plt.subplots(1,3, figsize=[18, 5])
            reco_countsx, reco_binsx, _x = axs[0].hist(df_filtered['reco_track'+pos+'x'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='ML Reco Track '+positions_labels[j])
            true_countsx, true_binsx, _x = axs[0].hist(df_filtered['true_track'+pos+'x'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='True Track '+positions_labels[j])
            axs[0].set_xlabel(positions_labels[j]+' Position [cm]', fontsize=12)
            axs[0].set_ylabel('Count [Tracks / 3 cm]', fontsize=12)
            axs[0].set_xlim(-75, 70)

            reco_countsy, reco_binsy, _y = axs[1].hist(df_filtered['reco_track'+pos+'y'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='ML Reco Track '+positions_labels[j])
            true_countsy, true_binsy, _y = axs[1].hist(df_filtered['true_track'+pos+'y'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='True Track '+positions_labels[j])
            axs[1].set_xlabel(positions_labels[j]+' Position [cm]', fontsize=12)
            axs[1].set_xlim(-65, 65)

            reco_countsz, reco_binsz, _z = axs[2].hist(df_filtered['reco_track'+pos+'z'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='ML Reco Track '+positions_labels[j])
            true_countsz, true_binsz, _z = axs[2].hist(df_filtered['true_track'+pos+'z'], alpha=0.5, \
                bins=100, range=(-150, 150), histtype='stepfilled', label='True Track '+positions_labels[j])
            axs[2].set_xlabel(positions_labels[j]+' Position [cm]', fontsize=12)
            axs[2].set_xlim(-75, 70)

            axs[1].legend(loc='upper right', fontsize=12)
            axs[0].set_title('       '+r'$\bf{X\ Coordinate}$', y=1.0, loc='left', pad=-14, fontsize=14)
            axs[1].set_title('       '+r'$\bf{Y\ Coordinate}$', y=1.0, loc='left', pad=-14, fontsize=14)
            axs[2].set_title('       '+r'$\bf{Z\ Coordinate}$', y=1.0, loc='left', pad=-14, fontsize=14)
            fig.suptitle('True vs. ML Reco Track '+positions_labels[j]+' Position for Reconstructed '+sample_label+' Sample', fontsize=16)

            plt.tight_layout()
            output.savefig()
            plt.close()

        # Plot Charged Track Multiplicity at Vertex
        fig = plt.figure(figsize=[8,6])
        axs = fig.add_subplot(111)
        reco_counts, reco_bins, _ = axs.hist(df_filtered['reco_ixn_charged_track_mult'], \
            alpha=0.5, bins=50, range=(0, 50), histtype='stepfilled', label='ML Reco Charged Track Multiplicity')
        true_counts, true_bins, _ = axs.hist(df_filtered['true_ixn_charged_track_mult'], \
            alpha=0.5, bins=50, range=(0, 50),histtype='stepfilled', label='True Charged Track Multiplicity')
        axs.set_xlabel('Charged Track Multiplicity at Vertex', fontsize=12)
        axs.set_ylabel('Count [Tracks / Track]', fontsize=12)
        axs.set_xlim(0, 20)
        axs.legend(loc='upper right', fontsize=12)
        axs.set_title('True vs. ML Reco Charged Track Multiplicity at Vertex for\n Reconstructed'+sample_label+'Sample', fontsize=15)
        plt.tight_layout()
        output.savefig()
        plt.close()





# Specific plots: plotting angles by length
def plot_angle_by_length(sample_name,angle_field,angle_label,df_length_grouped,\
                         bins_shown,max_bins,max_track_length):

    height = (bins_shown/max_bins)*45

    fig, axs = plt.subplots(bins_shown, sharex=True, figsize=[10, height])
    max_count = 50.
    for i, (name, group) in enumerate(df_length_grouped):
        if i>=bins_shown:
            break
        reco_counts, reco_bins, _ = axs[i].hist(abs(np.cos(group['reco_'+angle_field])), alpha=0.5, \
                                                bins=50, range=(0, 1), histtype='stepfilled',\
                                                label='ML Reco '+angle_label)
        true_counts, true_bins, _ = axs[i].hist(abs(np.cos(group['true_'+angle_field])), alpha=0.5, \
                                                bins=50, range=(0, 1), histtype='stepfilled',\
                                                label='True '+angle_label)
        if i==0:
            max_count = np.max(np.array([np.max(reco_counts), np.max(true_counts)]))
            print(max_count)
            axs[i].legend(loc='upper right', fontsize=12)
        axs[i].set_xlabel('| cos('+angle_label+') |', fontsize=12)
        axs[i].set_ylabel('Count ['+sample_name+' / 0.02]', fontsize=12)
        #axs[i].legend(loc='upper right')
        axs[i].set_xlim(0, 1)
        axs[i].set_yscale('log')
        axs[i].set_ylim(None,max_count+1000)
        length_min = max_track_length/max_bins*(name-1)
        length_max = max_track_length/max_bins*name
        axs[i].set_title('  '+rf'$\bf{{Reco\ Track\ Lengths:\ {length_min:.0f} - {length_max:.0f}\ cm}}$', y=1.0, loc='left', pad=-14, fontsize=14)
    fig.suptitle('True vs. ML Reco Absolute Value of Cosine of\n'+angle_label+' for Reconstructed '+sample_name+' Sample', fontsize=16)

    plt.tight_layout()
    plt.subplots_adjust(top=0.95, hspace=0.0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--ttree_file', default=None, type=str,\
                        help='''string corresponding to the path of the file containing a \
                                TTree with a sample of reconstructed track characteristics''')
    parser.add_argument('-p', '--protons_only', default=False, type=bool,\
                        help='''bool corresponding to whether or not plots will be for all charged \
                                tracks or just those associated with protons''')
    args = parser.parse_args()
    main(**vars(args))