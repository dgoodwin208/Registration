% Various visualiztion scripts which might be useful
% Useful only if the user puts a breakpoint after the calc_affine.m script
% runs to completion. As of Dec 2015, the breakpoint should be put in line 
% 218 of registerWithDescriptors.m
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

