# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:35:07 2024

@author: sf02
"""

# this script is to be run after make_coefficient_surface_mesh.sh and before make_coefficient_videos.sh.

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from nilearn import surface, plotting

analyses = ['dB_power','dB_theta_power','dB_alpha_power','dB_beta_power','dB_gamma_power','dB_high_gamma_power','voltage']

# set path
dirp='/group/mlr-lab/Saskia/ECoG_LASSO/'
sys.path.append(dirp)

# get pial surface mesh and curvature data
sys.path.append(dirp+'/scripts/fsaverage')
pial_left = dirp+'/scripts/fsaverage/pial_left.gii.gz'
curv_left = dirp+'/scripts/fsaverage/curv_left.gii.gz'
pial_right = dirp+'/scripts/fsaverage/pial_right.gii.gz'
curv_right = dirp+'/scripts/fsaverage/curv_right.gii.gz'

# set output directory
os.mkdir(dirp+'/derivatives/figures/coefficients/')
imagedir=dirp+'/derivatives/figures/coefficients/surface_plots/'
os.mkdir(imagedir)

# 17 patients, NOT 18, because we have no coordinates for patient 20.
# 14 left hemisphere and 3 right hemisphere 
Lsubcode = ['sub-01','sub-02','sub-03','sub-04','sub-05','sub-07','sub-09','sub-10','sub-11','sub-13','sub-14','sub-15','sub-21','sub-22']
Rsubcode = ['sub-06','sub-08','sub-12']

# for each analysis
for current_analysis in analyses:
    # for each timepoint
    for t in range(166):
        # initialise
        Ltotals = np.zeros(163842)
        Lanimatecoefs = np.zeros(163842)
        Rtotals = np.zeros(163842)
        Ranimatecoefs = np.zeros(163842)
        # for each left-hemisphere participant
        for s in Lsubcode:
            # Load coefficients
            Lfc = surface.load_surf_data(dirp+'/work/'+s+'/coefficients/surface/L_'+s+'_'+current_analysis+'_timepoint_{:03}_coefficients.func.gii'.format(t+1))
            tmp = Lfc != 0
            tmp = tmp.astype(int)
            Ltotals += tmp
            # update sum of NEGATIVE coefficients (N.B. Chris and Tim code animal as 1 and therefore count positive coefficients)
            tmp = Lfc < 0
            tmp = tmp.astype(int)
            Lanimatecoefs += tmp
        # for each right-hemisphere participant
        for s in Rsubcode:
            # do the same thing
            Rfc = surface.load_surf_data(dirp+'/work/'+s+'/coefficients/surface/R_'+s+'_'+current_analysis+'_timepoint_{:03}_coefficients.func.gii'.format(t+1))
            tmp = Rfc != 0
            tmp = tmp.astype(int)
            Rtotals += tmp
            tmp = Rfc < 0
            tmp = tmp.astype(int)
            Ranimatecoefs += tmp
        # set unselected vertices to nan
        Ltotals[Ltotals==0] = np.nan
        Rtotals[Rtotals==0] = np.nan
        # convert totals to proportion of participants in which that vertex is selected
        Lfinalselection = Ltotals/len(Lsubcode)
        Rfinalselection = Rtotals/len(Rsubcode)
        # plot selection
        # Left hemisphere viewed from underneath
        #fig = plotting.plot_surf_stat_map(
        #pial_left, Lfinalselection, hemi='left',view=[270,270],
        #colorbar=True, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
        #bg_on_data=True,
        #)
        #fig.savefig(imagedir+'/L_'+current_analysis+'_timepoint_{:03}_selection.jpeg'.format(t+1),dpi=300)
        #plt.close(fig)
        # Right hemisphere viewed from underneath
        #fig = plotting.plot_surf_stat_map(
        #pial_right, Rfinalselection, hemi='right',view=[270,270],
        #colorbar=True, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
        #bg_on_data=True,
        #)
        #fig.savefig(imagedir+'/R_'+current_analysis+'_timepoint_{:03}_selection.jpeg'.format(t+1),dpi=300)
        #plt.close(fig)
        
        # calculate coefficient directions
        Ldirections = Lanimatecoefs/Ltotals
        Rdirections = Ranimatecoefs/Rtotals
        # plot direction
        # Left hemisphere viewed from underneath
        fig = plotting.plot_surf_stat_map(
        pial_left, Ldirections, hemi='left',view=[270,270],
        colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
        bg_on_data=True,
        )
        fig.savefig(imagedir+'/L_'+current_analysis+'_timepoint_{:03}_coefficients.jpeg'.format(t+1),dpi=300)
        plt.close(fig)
        # Right hemisphere viewed from underneath
        fig = plotting.plot_surf_stat_map(
        pial_right, Rdirections, hemi='right',view=[270,270],
        colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
        bg_on_data=True,
        )
        fig.savefig(imagedir+'/R_'+current_analysis+'_timepoint_{:03}_coefficients.jpeg'.format(t+1),dpi=300)
        plt.close(fig)
        del(tmp,Lanimatecoefs,Ldirections,Lfc,Lfinalselection,Ltotals,Ranimatecoefs,Rdirections,Rfc,Rfinalselection,Rtotals)
    
