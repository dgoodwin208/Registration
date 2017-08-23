% Step 4-6: Calculating Correspondences 
% This is the code that calculates the keypoints and descriptors at 
% varying scale levels
%
% INPUTS: 
% moving_run: which expeirment do you want to warp accordingly? 
% OUTPUTS:
% no variables. All outputs saved to params.OUTPUTDIR
% 
% Author: Daniel Goodwin dgoodwin208@gmail.com 
% Date: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function registerWithDescriptors(moving_run)

    loadExperimentParams;

    params.MOVING_RUN = moving_run;

    disp(['RUNNING ON MOVING: ' num2str(params.MOVING_RUN) ', FIXED: ' num2str(params.FIXED_RUN) ', ' params.DATACHANNEL])

    filename = fullfile(params.INPUTDIR,sprintf('%sround%03d_%s.tif',...
                        params.SAMPLE_NAME,params.FIXED_RUN,params.DATACHANNEL));

    imgFixed_total = load3DTif(filename);


    filename = fullfile(params.INPUTDIR,sprintf('%sround%03d_%s.tif',...
                        params.SAMPLE_NAME,params.MOVING_RUN,params.DATACHANNEL));
    
    imgMoving_total = load3DTif(filename);

    %LOAD FILES WITH CROP INFORMATION, CROP LOADED FILES
    cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_cropbounds.mat',params.SAMPLE_NAME,params.FIXED_RUN));
    if exist(cropfilename,'file')==2
        load(cropfilename,'bounds'); bounds_fixed = floor(bounds); clear bounds;
        imgFixed_total = imgFixed_total(bounds_fixed(1):bounds_fixed(2),bounds_fixed(3):bounds_fixed(4),:);        
    end
    
    cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_cropbounds.mat',params.SAMPLE_NAME,params.MOVING_RUN));
    if exist(cropfilename,'file')==2
        load(cropfilename,'bounds'); bounds_moving = floor(bounds); clear bounds;
        imgMoving_total = imgMoving_total(bounds_moving(1):bounds_moving(2),bounds_moving(3):bounds_moving(4),:);
    end

    %------------------------------Load Descriptors -------------------------%
    %Load all descriptors for the MOVING channel
    descriptor_output_dir_moving = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s/',params.SAMPLE_NAME, ...
                                            params.MOVING_RUN,params.REGISTERCHANNEL));

    files = dir(fullfile(descriptor_output_dir_moving,'*.mat'));
    keys_moving_total = {}; keys_ctr=1;
    for file_idx= 1:length(files)
        filename = files(file_idx).name;
        
        %The data for each tile is keys, xmin, xmax, ymin, ymax
        data = load(fullfile(descriptor_output_dir_moving,filename));
        for idx=1:length(data.keys)
            try
            %copy all the keys into one large vector of cells
            keys_moving_total{keys_ctr} = data.keys{idx};
            keys_moving_total{keys_ctr}.x = data.keys{idx}.x + data.xmin-1;
            keys_moving_total{keys_ctr}.y = data.keys{idx}.y + data.ymin-1;
            
            keys_ctr = keys_ctr+ 1;
            catch
                fprintf('TEMPERROR in filename %s',filename);
            end
        end
    end
    
    %Load all descriptors for the FIXED channel
    descriptor_output_dir_fixed = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s/',params.SAMPLE_NAME, ...
                                            params.FIXED_RUN,params.REGISTERCHANNEL));
    %descriptor_output_dir_fixed = fullfile(params.OUTPUTDIR,sprintf('sample%dround%03d_%s/',params.SAMPLE_NUM,params.FIXED_RUN,params.REGISTERCHANNEL));

    files = dir(fullfile(descriptor_output_dir_fixed,'*.mat'));
    keys_fixed_total = {}; keys_ctr=1;

    for file_idx= 1:length(files)
        filename = files(file_idx).name;
        %All results from CalculateDescriptorsforTilesAtIndeices saves the
        %descriptors in the coords of the tile
        
        data = load(fullfile(descriptor_output_dir_fixed,filename));
        for idx=1:length(data.keys)
            %copy all the keys into one large vector of cells
            keys_fixed_total{keys_ctr} = data.keys{idx};        %#ok<*AGROW>
            keys_fixed_total{keys_ctr}.x = data.keys{idx}.x + data.xmin-1;
            keys_fixed_total{keys_ctr}.y = data.keys{idx}.y + data.ymin-1;
            
            
            keys_ctr = keys_ctr+ 1;
        end
    end

    %------------All descriptors are now loaded as keys_*_total -------------%



    %chop the image up into grid
    tile_upperleft_y_moving = floor(linspace(1,size(imgMoving_total,1),params.ROWS_TFORM+1));
    tile_upperleft_x_moving = floor(linspace(1,size(imgMoving_total,2),params.COLS_TFORM+1));

    %don't need to worry about padding because these tiles are close enough in
    %(x,y) origins
    tile_upperleft_y_fixed = floor(linspace(1,size(imgFixed_total,1),params.ROWS_TFORM+1));
    tile_upperleft_x_fixed = floor(linspace(1,size(imgFixed_total,2),params.COLS_TFORM+1));

    %loop over all the subsections desired for the piecewise affine, finding
    %all relevant keypoints then calculating the transform from there
    keyM_total = [];
    keyF_total = [];

    %Because it takes about 5-10 minutes to generate the global list of vetted
    %keys, after we generate them we now save them in the output_keys_filename
    %if it's aready been generated, we can skip directly to the TPS calculation
    output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_%sround%03d.mat',params.SAMPLE_NAME,params.MOVING_RUN));

    %If we need to run the robust model checking to identify correct
    %correspondences
    if ~exist(output_keys_filename,'file')
        
        for x_idx=1:params.COLS_TFORM
            for y_idx=1:params.ROWS_TFORM
                
                disp(['Running on row ' num2str(y_idx) ' and col ' num2str(x_idx) ]);
                
                %the moving code is defined the linspace layout of dimensions above
                ymin_moving = tile_upperleft_y_moving(y_idx);
                ymax_moving = tile_upperleft_y_moving(y_idx+1);
                xmin_moving = tile_upperleft_x_moving(x_idx);
                xmax_moving = tile_upperleft_x_moving(x_idx+1);
                
                tile_img_moving = imgMoving_total(ymin_moving:ymax_moving, xmin_moving:xmax_moving,:);
                
    %             outputfilename_moving = fullfile(descriptor_output_moving_dir, ...
    %                 [num2str(ymin_moving) '-' num2str(ymax_moving) '_' num2str(xmin_moving) '-' num2str(xmax_moving) '.mat']);
                
                %Before calculating any features, make sure the tile is not empty
                if checkIfTileEmpty(tile_img_moving,params.EMPTY_TILE_THRESHOLD)
                    disp('Sees the moving tile to be empty');
                    continue
                end
                
                %FindRelevant keys not only finds the total keypoints, but converts
                %those keypoints to the scope of the specific tile, not the global
                %position
                keys_moving = findRelevantKeys(keys_moving_total, ymin_moving, ymax_moving,xmin_moving,xmax_moving);
                
                
                %Loading the fixed tiles is detemined by some extra overlap between
                %the tiles (may not be necessary)
                tile_img_fixed_nopadding = imgFixed_total(tile_upperleft_y_fixed(y_idx):tile_upperleft_y_fixed(y_idx+1), ...
                    tile_upperleft_x_fixed(x_idx):tile_upperleft_x_fixed(x_idx+1), ...
                    :);
                tilesize_fixed = size(tile_img_fixed_nopadding);
                
                ymin_fixed = floor(max(tile_upperleft_y_fixed(y_idx)-(params.OVERLAP/2)*tilesize_fixed(1),1));
                ymax_fixed = floor(min(tile_upperleft_y_fixed(y_idx+1)+(params.OVERLAP/2)*tilesize_fixed(1),size(imgFixed_total,1)));
                xmin_fixed = floor(max(tile_upperleft_x_fixed(x_idx)-(params.OVERLAP/2)*tilesize_fixed(2),1));
                xmax_fixed = floor(min(tile_upperleft_x_fixed(x_idx+1)+(params.OVERLAP/2)*tilesize_fixed(2),size(imgFixed_total,2)));
                
                clear tile_img_fixed_nopadding;
                tile_img_fixed = imgFixed_total(ymin_fixed:ymax_fixed, xmin_fixed:xmax_fixed,:);
                
                if checkIfTileEmpty(tile_img_fixed,params.EMPTY_TILE_THRESHOLD)
                    disp('Sees the moving tile to be empty');
                    continue
                end
                
                %FindRelevant keys not only finds the total keypoints, but converts
                %those keypoints to the scope of the specific tile, not the global
                %position
                keys_fixed = findRelevantKeys(keys_fixed_total, ymin_fixed, ymax_fixed,xmin_fixed,xmax_fixed);
                
                
                disp(['Sees ' num2str(length(keys_fixed)) ' features for fixed and ' num2str(length(keys_moving)) ' features for moving.']);
                if length(keys_fixed)==0 || length(keys_moving)==0
                    disp('Empty set of descriptors. Skipping')
                    continue;
                end
                
                
                
                % ----------- SIFT MATCHING AND ROBUST MODEL SELECTION ----------%
                %
                DM = []; %M for moving
                LM = [];
                for i = 1:length(keys_moving)-1
                    DM(i,:) = keys_moving{i}.ivec;
                    LM(i,:) = [keys_moving{i}.y, keys_moving{i}.x, keys_moving{i}.z];
                end
                DF = []; %F for fixed
                LF = [];
                for i = 1:length(keys_fixed)-1
                    DF(i,:) = keys_fixed{i}.ivec;
                    LF(i,:) = [keys_fixed{i}.y, keys_fixed{i}.x, keys_fixed{i}.z];
                end
                correspondences = vl_ubcmatch(DM',DF');
                
                if length(correspondences)<20
                    disp(['We only see ' num2str(length(correspondences)) ' which is insufficient to calculate a reliable transform. Skipping']);
                    continue
                end
                
                try
                    calc_affine;
                    %calc_affine produces keyM and keyF, pairs of point correspondences
                    %from the robust model fitting. The math is done with local
                    %coordinates to the subvolume, so it needs to be adapted to global
                    %points

                    keyM_total = [keyM_total; keyM(:,1) + ymin_moving-1, keyM(:,2) + xmin_moving-1, keyM(:,3) ];
                    keyF_total = [keyF_total; keyF(:,1) + ymin_fixed-1, keyF(:,2) + xmin_fixed-1, keyF(:,3)];
                
                catch
                    disp(['Cannot process sub-image index: ' num2str(y_idx) ' ' num2str(x_idx)])
                    continue;
                end
                % ----------- END ---------- %
                
                
                
            end
        end
        
        save(output_keys_filename,'keyM_total','keyF_total');
    else %if we'va already calculated keyM_total and keyF_total, we can just load it
        disp('KeyM_total and KeyF_total already calculated. Loading.');
        load(output_keys_filename);
    end

    % Do we filter on maximum displacements in coordinate space?
    % For some data this might be useful
    if (params.MAXDISTANCE>-1)
        remove_indices = [];
        for match_idx = 1:size(keyF_total,1)
           if norm(keyF_total(match_idx,:)-keyM_total(match_idx,:))>params.MAXDISTANCE
            remove_indices = [remove_indices match_idx];
           end
        end
        keyF_total(remove_indices,:) = [];
        keyM_total(remove_indices,:) = [];
        clear remove_indices;
    end
    
    %Do a global affine transform on the data and keypoints before
    %doing the fine-resolution non-rigid warp
   
    %Because of the annoying switching between XY/YX conventions, 
    %we have to switch XY components for the affine calcs, then switch back
    keyM_total_switch = keyM_total(:,[2 1 3]);
    keyF_total_switch = keyF_total(:,[2 1 3]);
    affine_tform = findAffineModel(keyM_total_switch, keyF_total_switch);
    
    %Warp the keyM features into the new space
    rF = imref3d(size(imgFixed_total));
    %Key total_affine is now with the switched XY
    keyM_total_affine = [keyM_total_switch, ones(size(keyM_total_switch,1),1)]*affine_tform'; 
    %keyM_total is now switched 
    keyM_total=keyM_total_affine(:,1:3);
    %keyF_total = keyF_total_switch;
    %Remove any keypoints which are now outside the bounds of the image
    filtered_correspondence_indices = (keyM_total(:,1) <1 | keyM_total(:,2)<1 | keyM_total(:,3)<1 | ...
                                       keyM_total(:,1) > size(imgFixed_total,2) | ...
                                       keyM_total(:,2) > size(imgFixed_total,1) | ...
                                       keyM_total(:,3) > size(imgFixed_total,3) );
    fprintf('Losing %i features after affine warp\n',sum(filtered_correspondence_indices));
    keyM_total(filtered_correspondence_indices,:) = [];
    keyF_total(filtered_correspondence_indices,:) = [];
    
    %switch keyM back to the other format for the TPS calcs
    keyM_total = keyM_total(:,[2 1 3]);
        
    for c = 1:length(params.CHANNELS)
        %Load the data to be warped
        tic;
        data_channel = params.CHANNELS{c};
        fprintf('load 3D file for affine transform on %s channel\n',data_channel);
        filename = fullfile(params.INPUTDIR,sprintf('%sround%03d_%s.tif',params.SAMPLE_NAME,params.MOVING_RUN,data_channel));
        imgToWarp = load3DTif(filename);
        toc;
        
        output_affine_filename = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s_affine.tif',...
                        params.SAMPLE_NAME,params.MOVING_RUN,data_channel));
        imgMoving_total_affine = imwarp(imgToWarp,affine3d(affine_tform'),'OutputView',rF);
        save3DTif(imgMoving_total_affine,output_affine_filename);
    end
    
    disp('set up cluster and parpool')
    tic;
    cluster = parcluster('local_96workers');
    parpool(cluster,24)
    toc;

    output_TPS_filename = fullfile(params.OUTPUTDIR,sprintf('TPSMap_%sround%03d.mat',params.SAMPLE_NAME,params.MOVING_RUN));
    if exist(output_TPS_filename,'file')==0
%        [in1D_total,out1D_total] = TPS3DWarpWhole(keyM_total,keyF_total, ...
        
        [in1D_total,out1D_total] = TPS3DWarpWholeInParallel(keyM_total,keyF_total, ...
                            size(imgMoving_total), size(imgFixed_total));
        disp('save TPS file')
        tic;
        save(output_TPS_filename,'in1D_total','out1D_total','-v7.3');
        toc;
    else
        %load in1D_total and out1D_total
        disp('load TPS file')
        tic;
        load(output_TPS_filename);
        toc;
        %Experiments 7 and 8 may have been saved with zeros in the 1D vectors
        %so this removes it
        [ValidIdxs,I] = find(in1D_total>0);
        in1D_total = in1D_total(ValidIdxs);
        out1D_total = out1D_total(ValidIdxs);
    end


                        
    

    %Warp all three channels of the experiment once the index mapping has been
    %created
    for c = 1:length(params.CHANNELS)
        %Load the data to be warped
        disp('load 3D file to be warped')
        tic;
        data_channel = params.CHANNELS{c};
        filename = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s_affine.tif',params.SAMPLE_NAME,params.MOVING_RUN,data_channel));
        imgToWarp = load3DTif(filename);
        toc;
        
        %we loaded the bounds_moving data at the very beginning of this file
        if exist(cropfilename,'file')==2
            imgToWarp = imgToWarp(bounds_moving(1):bounds_moving(2),bounds_moving(3):bounds_moving(4),:);
        end
        [ outputImage_interp ] = TPS3DApply(in1D_total,out1D_total,imgToWarp,size(imgFixed_total));
        
        outputfile = fullfile(params.OUTPUTDIR,sprintf('%sround%03d_%s_registered.tif',params.SAMPLE_NAME,params.MOVING_RUN,data_channel));
        save3DTif(outputImage_interp,outputfile);
    end

    disp('delete parpool')
    tic;
    p = gcp('nocreate');
    delete(p);
    toc;







end


