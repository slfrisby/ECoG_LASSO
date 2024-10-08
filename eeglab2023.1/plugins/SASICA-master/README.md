This is SASICA, a tool to help you reject/select independent components based on various properties of these components.

This is a pre-release of the FieldTrip compatible version of SASICA (originally for EEGLAB). Please report any bug to max.chaumon@gmail.com

# Install:
In your favorite directory with `git clone -b feature/ft_compat https://github.com/dnacombo/SASICA.git`
or download the zip file [here](https://github.com/dnacombo/SASICA.git) and uncompress in your favorite directory.
Then in MATLAB: `addpath('your/favorite/directory')`


# Usage:

	[cfg] = ft_SASICA(cfg,comp,data)
	with inputs:
         - cfg: a structure with field
               layout: that will be passed to ft_prepare_layout
              as well as any of the following fields (see below for explanations):
               autocorr    detect components with low autocorrelation
               focalcomp   detect focal components in sensor space
               trialfoc    detect focal components in trial space
               SNR         detect components with low signal to noise
                           ratio across trials between two time windows.
               EOGcorr     detect components with high correlation with
                           vertical and horizontal EOG channels
               chancorr    detect components with high correlation with
                           any channel
               FASTER      use FASTER (Nolan et al. 2010) detection
                           methods.
               ADJUST      use ADJUST (Mongon et al. 2011) detection
                           methods
               opts        set various options: noplot, nocompute, FontSize


        - comp: the output of ft_componentanalysis
        - data: the output of ft_preprocessing
For more detailed information, see `doc eeg_SASICA`

For an example cfg structure, run `cfg = ft_SASICA('getdefs')`


# Available methods are:

- Autocorrelation: detects noisy components with weak 
    autocorrelation (muscle artifacts usually)
- Focal components: detects components that are too focal and 
    thus unlikely to correspond to neural 
    activity (bad channel or muscle usually).
- Focal trial activity: detects components with focal trial 
    activity, with same algorhithm as focal 
    components above. Results similar to trial 
    variability.
- Signal to noise ratio: detects components with weak signal 
    to noise ratio between arbitrary baseline 
    and interest time windows.
- EOG correlation: detects components whose time course 
    correlates with EOG channels.
- Bad channel correlation: detects components whose time course 
    correlates with any channel(s).
- ADJUST selection: use ADJUST routines to select components 
    (see Mognon, A., Jovicich, J., Bruzzone, 
    L., & Buiatti, M. (2011). ADJUST: An 
    automatic EEG artifact detector based on 
    the joint use of spatial and temporal 
    features. Psychophysiology, 48(2), 229-240. 
    doi:10.1111/j.1469-8986.2010.01061.x) 
- FASTER selection: use FASTER routines to select components
    (see Nolan, H., Whelan, R., & Reilly, R. B.
    (2010). FASTER: Fully Automated Statistical
    Thresholding for EEG artifact Rejection.
    Journal of Neuroscience Methods, 192(1),
    152-162. doi:16/j.jneumeth.2010.07.015)


 If you use this program in your research, please cite the following
 article:
 
> Chaumon M, Bishop DV, Busch NA. A Practical Guide to the Selection of
> Independent Components of the Electroencephalogram for Artifact
> Correction. Journal of neuroscience methods. 2015  


Copyright (C) 2019  Maximilien Chaumon

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Some of the measures used here are based on http://bishoptechbits.blogspot.com/2011/05/automated-removal-of-independent.html



