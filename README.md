# LiveInitiationDeconvolution
Analysis code for performing Bayesian deconvolution
------------------------------

The Matlab scripts and functions are provided to perform the deconvolution of initiation events and estimate bursting parameters as described in the manuscript. Test data is provided for hb NC14.

The code has been tested on MATLAB R2023b.

One can find in the repository the following piece of codes:
1. The main script to format the data, run the deconvolution and plot the output.
    - MainDeconvolution.m
2. Main functions that performs the deconvolution and the estimation of bursting parameters from transcriptional time traces.
    - DeconvSampling.m
    - ExtractIniPheno.m
3. Two utility functions.
    - makeKernel.m
    - plotEnveloppe.m
4. The input and output test data as zip files.
    - InDataSplit.zip
    - OutDeconvResSplit.zip

%   Copyright (c) 2024, Benjamin Zoller  
%   All rights reserved.  
%  
%   This source code is licensed under the MIT license found in the  
%   LICENSE file in the root directory of this source tree.