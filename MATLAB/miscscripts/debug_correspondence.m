% Various visualiztion scripts which might be useful
% Useful only if the user puts a breakpoint after the calc_affine.m script
% runs to completion. As of Dec 2015, the breakpoint should be put in line 
% 218 of registerWithDescriptors.m, before the correspondences have been
% calculated and saved as a file globalkeys_sample1round2.mat in
% params.OUTPUTDIR. If this file has been created, delete it for using this
% debug scripts

% Author: Daniel Goodwin dgoodwin208@gmail.com 
% Date: August 2015

%% View a point cloud of discovered keypoints
% This can be useful to make sure that the keypoints are reasonably
% distributed. Illumination artifacts from the microscope, for example, can
% cause clear areas of missing keypoints. Similarly, this will show
% subvolumes that were not processed in CalculateDscriptsForTiles.m

figure;
scatter3(LM(correspondences(1,1:num_params),1),LM(correspondences(1,1:num_params),2),LM(correspondences(1,1:num_params),3))
hold on;
scatter3(LM(correspondences(1,:),1),LM(correspondences(1,:),2),LM(correspondences(1,:),3),'r.')
hold off;
xlim([1 size(tile_img_fixed,2)]); xlabel('X coordinate');
ylim([1 size(tile_img_fixed,1)]); ylabel('Y coordinate');
zlim([1 size(tile_img_fixed,3)]); zlabel('Z coordinate');

%% View the correspondences between the images

for i = vote_inliers_idx(1:num_params)
    keyM = LM(correspondences(1,i),:);
    keyF = LF(correspondences(2,i),:);
    keyF_est = floor([keyM 1] * ransac_tform');
    
    disp([mat2str(keyM) ' ' mat2str(keyF_est(1:3)) ' ' mat2str(keyF)]); 
    %keyF = LF(correspondences(2,i),:);
    
    figure(1);
    subplot (1,2,1);
    hold on;
    img_slice_moving = tile_img_moving(:,:,keyM(3));
    imshow(img_slice_moving ,[min(min(img_slice_moving)),max(max(img_slice_moving))]);
    plot (keyM(1), keyM(2), 'b*');
    title(['Moving img. idx: ' num2str(correspondences(1,i)) ' value ' mat2str(keyM) ' ct ' num2str(inlier_stat_vec(i)) ])
    hold off;
    
    subplot (1,2,2);
    hold on;
    img_slice_fixed = tile_img_fixed(:,:,keyF(3));
    imshow (img_slice_fixed,[min(min(img_slice_fixed)),max(max(img_slice_fixed))]);
    plot (keyF(1), keyF(2), 'r*');
    title(['Fixed image. idx: ' num2str(correspondences(2,i)) ' value' mat2str(keyF) ' ct ' num2str(inlier_stat_vec(i)) ])
    hold off;
    pause
end


%% Make a plot of the best points by a voting method
common_inliers = find(inlier_stat_vec>100);

final_inliers = zeros(length(correspondences));
final_inliers(best_inliers) = 1;

[sorted_inliers, indices] = sort(inlier_stat_vec);
figure; plot(sorted_inliers);
hold on;
for i=1:length(indices) %looping over all indices
   %The sorted index is 
   sorted_pos = indices(i);
    
   if final_inliers(sorted_pos) ==1
    plot(i,inlier_stat_vec(indices(i)),'r*') 
   end
end
hold off;
title('Count of inlier frequency. Final Ransac model inliers highlighted')
xlabel('Sorted indices')
ylabel('Inlier count after RANSAC')

%% Create an figure that shows the different subregions per image
%chop the image up into grid
loadExperimentParams;
loadandCrop;

tile_upperleft_y_moving = floor(linspace(1,size(imgMoving_total,1),params.ROWS_TFORM+1));
tile_upperleft_x_moving = floor(linspace(1,size(imgMoving_total,2),params.COLS_TFORM+1));

%don't need to worry about padding because these tiles are close enough in
%(x,y) origins
tile_upperleft_y_fixed = floor(linspace(1,size(imgFixed_total,1),params.ROWS_TFORM+1));
tile_upperleft_x_fixed = floor(linspace(1,size(imgFixed_total,2),params.COLS_TFORM+1));

%% temp
figure(1)
subplot(params.COLS_TFORM,params.ROWS_TFORM,1)
figure(2)
subplot(params.COLS_TFORM,params.ROWS_TFORM,1)

for x_idx=1:params.COLS_TFORM
            for y_idx=1:params.ROWS_TFORM
                
                %the moving code is defined the linspace layout of dimensions above
                ymin_moving = tile_upperleft_y_moving(y_idx);
                ymax_moving = tile_upperleft_y_moving(y_idx+1);
                xmin_moving = tile_upperleft_x_moving(x_idx);
                xmax_moving = tile_upperleft_x_moving(x_idx+1);
                
                tile_img_moving = imgMoving_total(ymin_moving:ymax_moving, xmin_moving:xmax_moving,:);
                
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
                
                %now make the plots:
                figure(1) %fixed
                subplot(params.COLS_TFORM,params.ROWS_TFORM,(y_idx-1)*params.COLS_TFORM+x_idx)
                imagesc(tile_img_fixed(:,:,size(tile_img_fixed,3)/2)); axis off
                figure(2) %moving
                subplot(params.COLS_TFORM,params.ROWS_TFORM,(y_idx-1)*params.COLS_TFORM+x_idx)
                imagesc(tile_img_moving(:,:,size(tile_img_moving,3)/2)); axis off;
                
            end
end
