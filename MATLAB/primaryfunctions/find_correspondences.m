%Convert the 3D sift features into a usable array
function correspondences = find_correspondences(keys_moving, keys_fixed)

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

    correspondences = vl_ubcmatch(DM',DF',1.5); %1.5 is default

end
% %% Print and plot the results for sanity
% 
% template = '(%i,%i,%i) and (%i,%i,%i) \n';
% 
% for i = 1:length(correspondences)
%     locM = LM(correspondences(1,i),:);
%     locF = LF(correspondences(2,i),:);
%     fprintf(template,locM(1),locM(2),locM(3),locF(1),locF(2),locF(3) );
% end
% 
% disp(['For ' num2str(length(DM)) ' keypoints, found ' num2str(length(correspondences)) ' matches'])
% 
% % Pic a z-slice for the first image, then find all the matches to the second
% subplot (1,2,1);
% z1 = 226;
% 
% keypointsM = [];
% keypointsF = [];
% 
% zF_list = [];
% for idx = 1:length(correspondences)
%     locM_idx = correspondences(1,idx);
%     locF_idx = correspondences(2,idx);
%     
%     if LM(locM_idx,3)==z1
%         keypointsM = [keypoints1; L1(loc1_idx,1:2)];
%         keypointsF = [keypoints2; L2(loc2_idx,1:2)];
%         zF_list = [zF_list; LF(locF_idx,3)];
%     end
% end
% 
% % keypointsM
% % keypointsF
% % zF_list 
% 
% zF = floor(mean(zF_list));
% 
% figure(1)
% subplot (1,2,1);
% hold on;
% img_slice_moving = test_img_moving(:,:,zF);
% imshow(img_slice_moving ,[min(min(img_slice_moving)),max(max(img_slice_moving))]);
% plot (keypointsM(:,2), keypointsM(:,1), 'b*');
% hold off;
% 
% subplot (1,2,2);
% hold on;
% img_slice_fixed = test_img_fixed(:,:,zF);
% imshow (test_img_fixed(:,:,zF),[min(min(img_slice_fixed)),max(max(img_slice_fixed))]);
% plot (keypointsF(:,2), keypointsF(:,1), 'r*');
% hold off;