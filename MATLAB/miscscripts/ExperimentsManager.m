% EXPERIMENTSMANAGER Goes through the input and output directories 
%   Detailed explanation goes here
% 
% First, get the experiment numbers

function [experiments,crops,keypts,registrations,cell_detections] = ExperimentsManager()
loadExperimentParams;

% filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
%     params.SAMPLE_NAME,params.MOVING_RUN,params.REGISTERCHANNEL));

% files = dir(fullfile(params.INPUTDIR,sprintf('%sround%d*.tif',params.SAMPLE_NAME,sample_num)));
files = dir(fullfile(params.INPUTDIR,sprintf('%sround%d*.tif',...
    params.SAMPLE_NAME,sample_num)));

experiments = [];
for file_idx = 1:length(files)
    filename = files(file_idx).name;
    underscore = strsplit(filename,'_'); 
    rounds = strsplit(underscore{2},'round');
    experiments(file_idx) = str2num(rounds{2});
end

experiments = unique(experiments);

crops = [];
keypts = [];
registrations = [];
cell_detections = [];

for exp_idx= 1:length(experiments)
    filename = fullfile(params.OUTPUTDIR, sprintf('%s%dround%d_cropbounds.mat', params.SAMPLE_NAME,sample_num, experiments(exp_idx)));
    crops(exp_idx) = exist(filename,'file') >0;
end

for exp_idx= 1:length(experiments)
    
    %filename = fullfile(params.OUTPUTDIR, sprintf('globalkeys_sample%dround%d.mat', sample_num, experiments(exp_idx)));
    filename = fullfile(params.OUTPUTDIR, sprintf('sample%dround%d_Lectin', sample_num, experiments(exp_idx)));
    
    keypts(exp_idx) = exist(filename,'dir') >0;
end

for exp_idx= 1:length(experiments)
    filenameDAPI = fullfile(params.OUTPUTDIR, sprintf('TPSsample%dround%d_DAPI', sample_num, experiments(exp_idx)));
    filenameData = fullfile(params.OUTPUTDIR, sprintf('TPSsample%dround%d_Data', sample_num, experiments(exp_idx)));
    filenameLectin = fullfile(params.OUTPUTDIR, sprintf('TPSsample%dround%d_Lectin', sample_num, experiments(exp_idx)));
    registrations(exp_idx) = (exist(filenameDAPI,'dir') >0) & (exist(filenameData,'dir') >0) & (exist(filenameLectin,'dir') >0) ;
end

for exp_idx= 1:length(experiments)
    filename = fullfile(params.OUTPUTDIR, sprintf('centroids_sample%dround%d_DAPI.mat', sample_num, experiments(exp_idx)));
    cell_detections(exp_idx) = exist(filename,'file') >0;
end


%Pretty print the results

disp(sprintf('ExpNum\tCropped\tKeypts\tRegist\tCell Ct'));

for i= 1:length(experiments)
    disp(sprintf('%d\t%d\t%d\t%d\t%d',experiments(i),crops(i),keypts(i),registrations(i),cell_detections(i)));
end
