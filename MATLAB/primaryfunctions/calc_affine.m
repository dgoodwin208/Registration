% This is ran inside the registerWithDescriptors function
% left as a script just to keep the code separate.

% Author: Daniel Goodwin, Aug 2015 dgoodwin208@gmail.com


%Create new A, calcuate the transformed values, get the distance from the
%target points as defined by the correspondences
keypointsM = [];
keypointsF = [];

for idx = 1:length(correspondences)
    locM_idx = correspondences(1,idx);
    locF_idx = correspondences(2,idx);
    
    keypointsM = [keypointsM; LM(locM_idx,:)];
    keypointsF = [keypointsF; LF(locF_idx,:)];   
end

%RANSAC params
threshold = 3.; %How close does it have to be??
best_score = -1;
best_inliers = [];
best_tform = [];
p = 0.99;         % Desired probability of choosing at least one sample
                  % free from outliers 
                  
inlier_stat_vec = zeros(1,length(correspondences));
trial_iter_count = 0;
N=10;  %initting N, will be set first iteration    
max_trial_limit = 2000000;
min_trial_limit = 10000;
%input points
num_inputs = length(correspondences);
num_samples = 6; %4 or any other number used to calculate the affine tform

while trial_iter_count < N || trial_iter_count<min_trial_limit
    
    permuted_indices = randperm(num_inputs);

    %randomly choose num_inputs from keyM and keyF
    random_indices = permuted_indices(1:num_samples);
    keyM = keypointsM(random_indices,:);
    keyF = keypointsF(random_indices,:);

    %Fit the model for points in Moving image to points in Fixed image 
    affine_tform = findAffineModel(keyM, keyF);

    inliers = [];
    
    for i = 1:num_inputs
        %The actual keypoint in the fixed image
        key_act = keypointsF(i,:);

        %apply model to the moving point to the get the point in fixed
        %coords
        homogeneous_res = [keypointsM(i,:) 1]*affine_tform';    
        
        key_est = homogeneous_res(1:3);    
        %disp([mat2str(key_act) ' ' mat2str(key_est)])
        
        dist = sqrt((key_act(1)-key_est(1))^2 + (key_act(2)-key_est(2))^2 + (key_act(3)-key_est(3))^2 );
        if dist <= threshold
           inliers = [inliers; i]; 
        end
    end

    ninliers = length(inliers);

    if ninliers > best_score
        best_score = ninliers;
        best_tform = affine_tform;
        best_inliers = inliers;

        % Update estimate of N, the number of trials to ensure we pick,
        % with probability p, a data set with no outliers.
        fracinliers =  ninliers/num_inputs;
        pNoOutliers = 1 -  fracinliers^4; %to the 4th bc minimum num points to make model
        pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
        pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
        N = log(1-p)/log(pNoOutliers);

        trial_iter_count = 0;

    end
    
    %update the vector to keep track of which correspondences are legit
    for i = 1:length(inliers)
       %If the index was part of the model calculation, then skip it
       %otherwise, if the point was an inlier, then add the count for that
       %time where the point is a legitimate inlier
       if ~any(inliers(i)==permuted_indices(1:num_samples))
           inlier_stat_vec(inliers(i)) = inlier_stat_vec(inliers(i)) +1;
       end
    end
    
    
    trial_iter_count = trial_iter_count +1;
    
    if trial_iter_count > max_trial_limit
        disp('maxing out of iterations')
        break
    end

end

disp([' best_score ' num2str(best_score) ' / ' num2str(num_inputs) ])

% To conclude, recalculate with the best point
% Calculate the transform with the inliers:
keyM = keypointsM(best_inliers,:);
keyF = keypointsF(best_inliers,:);

% Fit the model for points in Moving image to points in Fixed image 
affine_tform_tot = findAffineModel(keyM, keyF);

% disp(['Key Actual | Key Estimated | Original Keypoint Moving'])
% for i = 1:length(best_inliers)
%     % The actual keypoint in the fixed image
%     key_act = keypointsF(best_inliers(i),:);
% 
%     % apply model to the moving point to the get the point in fixed
%     % coords
%     homogeneous_res = [keypointsM(best_inliers(i),:) 1]*affine_tform_tot';    
% 
%     key_est = homogeneous_res(1:3);    
%     disp([mat2str(key_act) ' ' mat2str(floor(key_est)) ' ' mat2str(keypointsM(best_inliers(i),:)) ])
% end

ransac_tform = affine_tform_tot;



% Loop over all correspondences to visualize what a high voted, low voted correspondence looks like
%  What happens if we create a model from the most common inliers?
[x, vote_inliers_idx] = sort(inlier_stat_vec,'descend');
%Set the number of parameters to be all of those that are at least 10% of
%the max voted entry (done to avoid incorrect matches)
num_params = length(find(inlier_stat_vec>.1*max(inlier_stat_vec))); %floor(.25*length(vote_inliers_idx)); 
keyM = LM(correspondences(1,vote_inliers_idx(1:num_params)),:);
keyF = LF(correspondences(2,vote_inliers_idx(1:num_params)),:);



return




