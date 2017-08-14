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

function calculateDescriptors(run_num,varargin)
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
filename = fullfile(params.INPUTDIR,sprintf('%sround%03d_%s.tif',...
    params.SAMPLE_NAME,run_num,params.REGISTERCHANNEL));
img = load3DTif(filename);


cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_cropbounds.mat',params.SAMPLE_NAME,run_num));
if exist(cropfilename,'file')==2
    load(cropfilename,'bounds');
    img = img(bounds(1):bounds(2),bounds(3):bounds(4),:);
else
    fprintf('NOTE: No _cropbounds file found, so using the entire image. \nTo non-destructively crop your image,use crop_sample.m \n');
end

%chop the image up into grid
tile_upperleft_y = floor(linspace(1,size(img,1),params.ROWS_DESC+1));
tile_upperleft_x = floor(linspace(1,size(img,2),params.COLS_DESC+1));

% Commenting out this optimization from Atsushi, only because I don't
% understand how/why this beneficial

% img_cache = containers.Map('KeyType','int32','ValueType','any');
%
% if length(target_indices) < params.ROWS_DESC*params.COLS_DESC*0.5 % magic number
%     for x_idx=1:params.COLS_DESC
%         for y_idx=1:params.ROWS_DESC
%             tile_counter = (x_idx-1)*params.ROWS_DESC+y_idx;
%
%             % only run kypts+descriptors for specified indices
%             if ~ismember(tile_counter,target_indices)
%                 continue
%             end
%
%             % get region, indexing column-wise
%             ymin = tile_upperleft_y(y_idx);
%             ymax = tile_upperleft_y(y_idx+1);
%             xmin = tile_upperleft_x(x_idx);
%             xmax = tile_upperleft_x(x_idx+1);
%
%             % create overlap region for calcuating features
%             % will remove all points in the overlap region after calculation
%             % but this avoids edge effects on any boundaries of
%             ymin_overlap = floor(max(tile_upperleft_y(y_idx)-(params.OVERLAP/2)*(ymax-ymin),1));
%             ymax_overlap = floor(min(tile_upperleft_y(y_idx+1)+(params.OVERLAP/2)*(ymax-ymin),size(img,1)));
%             xmin_overlap = floor(max(tile_upperleft_x(x_idx)-(params.OVERLAP/2)*(xmax-xmin),1));
%             xmax_overlap = floor(min(tile_upperleft_x(x_idx+1)+(params.OVERLAP/2)*(xmax-xmin),size(img,2)));
%
%             %Calculate the features on the larger (overlapping) regions
%             img_cache(tile_counter) = img(ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap,:);
%
%         end
%     end
%
%     clearvars img;
% end


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
        descriptor_output_dir = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s/',params.SAMPLE_NAME,run_num,params.REGISTERCHANNEL));
        if exist(descriptor_output_dir,'dir')==0
            mkdir(descriptor_output_dir);
        end
        
        % get region, indexing column-wise
        ymin = tile_upperleft_y(y_idx);
        ymax = tile_upperleft_y(y_idx+1);
        xmin = tile_upperleft_x(x_idx);
        xmax = tile_upperleft_x(x_idx+1);
        
        %Calculate the features on the larger (overlapping) regions
        %         if length(img_cache) > 0
        %             tile_img = img_cache(tile_counter);
        %         else
        % create overlap region for calcuating features
        % will remove all points in the overlap region after calculation
        % but this avoids edge effects on any boundaries of
        ymin_overlap = floor(max(tile_upperleft_y(y_idx)-(params.OVERLAP/2)*(ymax-ymin),1));
        ymax_overlap = floor(min(tile_upperleft_y(y_idx+1)+(params.OVERLAP/2)*(ymax-ymin),size(img,1)));
        xmin_overlap = floor(max(tile_upperleft_x(x_idx)-(params.OVERLAP/2)*(xmax-xmin),1));
        xmax_overlap = floor(min(tile_upperleft_x(x_idx+1)+(params.OVERLAP/2)*(xmax-xmin),size(img,2)));
        
        tile_img = img(ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap,:);
        %         end
        
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
            %keys = SWITCH_tile_processing(tile_img);
            keys = SWITCH_tile_processingInParallel(tile_img);
        end
        
        %There is a different terminology for the x,y coordinates that
        %needs to be noted. The 3DSIFT code treats x as the 0th
        %index (ie, the vertical dimension) whereas I am used to saying y
        %is the vertical direction.  This next bit of code both switches
        %the x and y back to my convention
        
        %Any coordinates less than these can be discarded
        margin_xmin = xmin - xmin_overlap;
        margin_ymin = ymin - ymin_overlap;
        %These are the coords above the normal xmax:xmin, ymin:ymax bounds
        margin_xmax = size(tile_img,2)-(xmax_overlap-xmax);
        margin_ymax = size(tile_img,1)-(ymax_overlap-ymax);
        
        indices_to_remove = [];
        for key_idx=1:length(keys)
            temp_key = keys{key_idx};
            temp_key.y = keys{key_idx}.x;
            temp_key.x = keys{key_idx}.y;
            
            %If the point is taken from the overlapping region, we can
            %ignore. This is mainly just to avoid edge effects in
            %calculating keypoints and descriptor
            if temp_key.y <= margin_ymin || temp_key.x <= margin_xmin || ...
                    temp_key.y > margin_ymax || temp_key.x > margin_xmax
                
                indices_to_remove = [indices_to_remove key_idx];
            else
                %set the new key coords:
                %Keep all the points relative to xmin and xmax by
                %subtracting the margin_xmin and margin_xmax
                %this is because registerWithDescriptors uses xmin and xmax
                keys{key_idx}.y = temp_key.y - margin_ymin;
                keys{key_idx}.x = temp_key.x - margin_xmin;
                
                if keys{key_idx}.y>(ymax-ymin+1)|| keys{key_idx}.x>(xmax-xmin+1)
                   error('Something was not right with the indexing in the tile');
                end
            end
        end
        
        %Remove any keys that were in the overlapping region.
        %         keys{indices_to_remove} = [];
        %         fprintf('Removed %i keypoints from the overlapping region\n',length(indices_to_remove));
        final_indices = ones(length(keys),1);
        final_indices(indices_to_remove)=0;
        keys = keys(logical(final_indices));
        fprintf('Removed %i/%i keypoints from the excess overlapping region\n',length(indices_to_remove),length(keys));
        
        
        save(outputfilename,'keys','ymin','xmin','ymax','xmax', 'params','run_num',...
            'ymin_overlap','ymax_overlap', 'xmin_overlap','xmax_overlap');
        
        clear keys;
        
    end
end
