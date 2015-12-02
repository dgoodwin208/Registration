% Step 2-3: Keypoints and Descriptors
% Calculate the SIFT Descriptors for a subset of subvolumes of the experiments
% as dictated in the loadExperimentParams.m script.
%
% INPUTS:
% sample_num is the index of the tissue sample to register
% run_num is the index of the experiment for the specified sample
% start_idx and end_idx specify the indices of subvolumes to calculate
%        keypoints and descriptors
%
% OUTPUTS:
% The output is a file that looks like <ymin>_<ymax>_<xmin>_<xmax>.m describing
% the pixel indexs of the subregion (example: 1150-1724_1-546.mat). These files
% contain the keypoints and descriptors in the keys cell array.
%
% These files are then later loaded by the registerWithDescriptors.m file
% Author: Daniel Goodwin dgoodwin208@gmail.com
% Date: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalculateDescriptorsForTileAtIndices(run_num,varargin)
% CALCULATEDESCRIPTORFORTILEATINDICES  Calculates keypoints and descriptors
%   CALCULATEDESCRIPTORFORTILEATINDICES(run_num) calculates across all
%   tiles (subvolumes)
%   CALCULATEDESCRIPTORFORTILEATINDICES(run_num, start_idx, end_idx)
%   calculates across just the tiles specified by the indices of start_idx
%   and end_idx. The partitioning details are established int he
%   loadExperimentParams.m
%
%   See also LOADEXPERIMENTPARAMS, PLUS.

%Load all the parameters per file
loadExperimentParams;

if size(varargin,2)==1
    error('ERROR: You must specify both a start and an end index');
elseif size(varargin,2)==2
    start_idx = varargin{1}; end_idx = varargin{2};
else
    start_idx = 1; end_idx = params.COLS_DESC*params.ROWS_DESC;
end
    
    
target_indices = start_idx:end_idx;



%Loading the tif file associated with the reference channel (ie,
%Lectin) for the image specified by run_num
filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
    params.SAMPLE_NAME,run_num,params.REGISTERCHANNEL));
img = load3DTif(filename);


cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%d_cropbounds.mat',params.SAMPLE_NAME,run_num));
if exist(cropfilename,'file')==2
    load(cropfilename,'bounds');
    img = img(bounds(1):bounds(2),bounds(3):bounds(4),:);
else
    fprintf('NOTE: No _cropbounds file found, so using the entire image. \nTo non-destructively crop your image,use crop_sample.m \n');
end

%chop the image up into grid
tile_upperleft_y = floor(linspace(1,size(img,1),params.ROWS_DESC+1));
tile_upperleft_x = floor(linspace(1,size(img,2),params.COLS_DESC+1));

tile_counter = 0; %create a manual counter to be used in the data partitioning

for x_idx=1:params.COLS_DESC
    for y_idx=1:params.ROWS_DESC
        tile_counter = tile_counter+1;
        
        % only run kypts+descriptors for specified indices
        if ~ismember(tile_counter,target_indices)
            continue
        end
        
        disp(['Running on row ' num2str(y_idx) ' and col ' num2str(x_idx) ]);
        
        %Make sure the folders for the descriptor outputs exist:
        descriptor_output_dir = fullfile(params.OUTPUTDIR,sprintf('%sround%d_%s/',params.SAMPLE_NAME,run_num,params.REGISTERCHANNEL));
        if exist(descriptor_output_dir,'dir')==0
            mkdir(descriptor_output_dir);
        end
        
        % get region, indexing column-wise
        ymin = tile_upperleft_y(y_idx);
        ymax = tile_upperleft_y(y_idx+1);
        xmin = tile_upperleft_x(x_idx);
        xmax = tile_upperleft_x(x_idx+1);
        tile_img = img(ymin:ymax, xmin:xmax,:);
        
        outputfilename = fullfile(descriptor_output_dir, ...
            [num2str(ymin) '-' num2str(ymax) '_' num2str(xmin) '-' num2str(xmax) '.mat']);
        
        %Before calculating any features, make sure the tile is not empty
        if checkIfTileEmpty(tile_img,params.EMPTY_TILE_THRESHOLD)
            disp('This subregion is empty. Skipping');
            continue
        end
        
        
        if exist(outputfilename,'file')>0 %Make sure that the descriptors have been calculated!
            continue;
        else
            keys = SWITCH_tile_processing(tile_img);
        end
        
        %There is a different terminology for the x,y coordinates that
        %needs to be noted. The 3DSIFT code treats x as the 0th
        %index (ie, the vertical dimension) whereas I am used to saying y
        %is the vertical direction.  This next bit of code both switches
        %the x and y back to my convention while also putting
        %the discovered keypoints back into global coodinates
        
        for key_idx=1:length(keys)
            temp_key = keys{key_idx};
            temp_key.y = keys{key_idx}.x;
            temp_key.x = keys{key_idx}.y;
            
            %set the new key coords:
            keys{key_idx}.y = temp_key.y;
            keys{key_idx}.x = temp_key.x;
            
        end
        
        save(outputfilename,'keys','ymin','xmin','ymax','xmax', 'params','run_num');
        
        clear keys;
        
    end
end