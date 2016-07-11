
params.SAMPLE_NAME = 'ROI3';
params.FIXED_RUN = 1;
params.MOVING_RUN = 2;

params.DATACHANNEL = '561';
params.REGISTERCHANNEL ='488';

%In the Murray et al 2015 this was {'Lectin', 'DAPI', 'Data}
params.CHANNELS = {'488','561','640'}; 

params.INPUTDIR = '../../input/'; 
params.OUTPUTDIR = '../../output/';

params.OVERLAP = .1; %used to be 10%

%how many subsections to calculate the descriptors?
params.ROWS_DESC = 3;
params.COLS_DESC = 3;
%over how many subsections to calculate the piecewise affine transforms?
params.ROWS_TFORM = 1;
params.COLS_TFORM = 1;

%if more than this percentage of pixels doesn't change in a tile
%consider it an empty tile and skip it
params.EMPTY_TILE_THRESHOLD = .85;

% SCALE_PYRAMID: Create a set of blurs to create a "scale pyramid", or in 
% non-computervision speak:
% Blurring by a progressive set of gaussians to adjust for scale 
% differences between SWITCH Experiments.    
% Using the documentation in fspecial3.m, the standard deviations of 
% the gaussian is defined as SIZE/2/2.354 so that FWHM equals half filter 
% size (http://en.wikipedia.org/wiki/FWHM).The blur_size values are chosen 
% with assumptions of the width of the Lectin vessels 
% (5-20 pixels at 2um/pix res) observed and assumption of minor scale 
% disparities (ie <20%)
params.SCALE_PYRAMID = 5:5:20;
