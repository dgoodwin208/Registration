%% Load features found for the FIXED and MOVING images

loadExperimentParams;
loadandCrop;
%------------------------------Load Descriptors -------------------------%
%Load all descriptors for the MOVING channel
descriptor_output_dir_moving = fullfile(params.OUTPUTDIR,sprintf('%sround%d_%s/',params.SAMPLE_NAME, ...
    params.MOVING_RUN,params.REGISTERCHANNEL));

files = dir(fullfile(descriptor_output_dir_moving,'*.mat'));
keys_moving_total = {}; keys_ctr=1;
for file_idx= 1:length(files)
    filename = files(file_idx).name;
    
    %The data for each tile is keys, xmin, xmax, ymin, ymax
    data = load(fullfile(descriptor_output_dir_moving,filename));
    for idx=1:length(data.keys)
        
        %copy all the keys into one large vector of cells
        keys_moving_total{keys_ctr} = data.keys{idx};
        keys_moving_total{keys_ctr}.x = data.keys{idx}.x + data.xmin-1;
        keys_moving_total{keys_ctr}.y = data.keys{idx}.y + data.ymin-1;
        
        keys_ctr = keys_ctr+ 1;
    end
end

%Load all descriptors for the FIXED channel
descriptor_output_dir_fixed = fullfile(params.OUTPUTDIR,sprintf('%sround%d_%s/',params.SAMPLE_NAME, ...
    params.FIXED_RUN,params.REGISTERCHANNEL));
%descriptor_output_dir_fixed = fullfile(params.OUTPUTDIR,sprintf('sample%dround%d_%s/',params.SAMPLE_NUM,params.FIXED_RUN,params.REGISTERCHANNEL));

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


%% Check out the features for a specific image
keys = keys_moving_total; %create a new variable for fixed/moving

%Pre-package the locations of the images by z
z_location_lookup = {};
for z=1:size(imgMoving_total,3) % change the looping param to depending 
    z
    relevants = [];ctr = 1;
    for k = 1:length(keys)
        if keys{k}.z<= z+2 && keys{k}.z >= z-2
            relevants(ctr) = k;
            ctr = ctr+ 1;
        end
    end
    z_location_lookup{z}=relevants;
end
%% 
figure;
for z=50:size(imgMoving_total,3) % change the looping param to depending 
                                % on the image we're exploring
    imagesc(imgMoving_total(:,:,z))
    hold on;
    relevant_keys = z_location_lookup{z};
    for k = relevant_keys
        plot(keys{k}.x,keys{k}.y,'g*')
    end
    
    hold off;
    title (sprintf('%d relevant keypoints at z=%d',length(relevant_keys),z))
    pause
end

%% Make a 3D plot of all features found
keys3D = zeros(length(keys_moving_total),3);
for k = 1:length(keys_moving_total)
   keys3D(k,:) = [keys_moving_total{k}.x,keys_moving_total{k}.y, keys_moving_total{k}.z];
end
plot3(keys3D(:,1),keys3D(:,2),keys3D(:,3),'.')


%% Exploring affine transforms on subregions
% place a breakpoint on line 130 in calc_affine to generate this
figure;
subplot(2,2,1)
imagesc(tile_img_fixed(:,:,size(tile_img_fixed,3)/2))
%imagesc(max(tile_img_fixed,[],3))
title('Fixed image')
subplot(2,2,2)
imagesc(tile_img_moving(:,:,size(tile_img_moving,3)/2))
title('Moving image')

tformObj = affine3d(affine_tform_tot');
Rcb = imref3d(size(tile_img_fixed));
[img_warped,rb] = imwarp(tile_img_moving,tformObj,'OutputView',Rcb, 'Interp', 'cubic');


img_out = zeros(size(tile_img_fixed,1),size(tile_img_fixed,2),3);

keypts_mv = LM(correspondences(1,vote_inliers_idx(1:num_params)),:);
keypts_fx = LF(correspondences(2,vote_inliers_idx(1:num_params)),:);

for idx = 1:length(keypts_fx)
    
    %I think we want to plot LM 2,1
    kF = keypts_fx(idx,:);
    kM = keypts_mv(idx,:);

    subplot(2,2,1)
    imagesc(tile_img_fixed(:,:,kF(3)));
    hold on;
    plot(kF(2),kF(1),'r*');
    hold off;
    title('Fixed image')
    
    subplot(2,2,2)
    imagesc(tile_img_moving(:,:,kM(3)));
    hold on;
    plot(kM(2),kM(1),'r*');
    hold off;
    title('Moving image')

    subplot(2,2,3)
    z = kF(3);
    img_out(:,:,1) = squeeze(img_warped(:,:,z));
    img_out(:,:,1) = img_out(:,:,1)/max(max(img_out(:,:,1)));
    img_out(:,:,2) = squeeze(tile_img_fixed(:,:,z));
    img_out(:,:,2) = img_out(:,:,2)/max(max(img_out(:,:,2)));

    imagesc(img_out)

    pause
end

%% Load matched features
output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_%sround%d.mat',params.SAMPLE_NAME,params.MOVING_RUN));
load(output_keys_filename);

figure;
scatter3(keyF_total(:,1),keyF_total(:,2),keyF_total(:,3))
hold on;
scatter3(keyM_total(:,1),keyM_total(:,2),keyM_total(:,3))
for i=1:size(keyF_total,1)
    
line3 = [keyF_total(i,:); keyM_total(i,:)];
plot3(line3(:,1),line3(:,2),line3(:,3),'k');
end
hold off;    