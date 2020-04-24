function derivative_boost_ycgosu(dats, SPMinfo, outputdir, varargin)
% apply derivative boost or estimate time to peak
% Basis functions should be orthogonalized before running this function!!
% 
%
% inputs:
%   dats: beta values in .nii ([number of basis function to boost X directory of each files])
%
%   SPMdirectory: directory of SPM.mat / super important because accor
%   
%   ouputdir: specify your output directory.
%
%   'int_cols', [int_cols in numerical indx]     
%   : columns of interest should be entered. This should be regressors of
%   betas in 'dats'.
%   
%   
% 