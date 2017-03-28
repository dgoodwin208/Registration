% Step 4: Keypoints and Descriptors
% This is the code that calculates the keypoints and descriptors at 
% varying scale levels
%
% INPUTS: 
% img is the image volume
% OUTPUTS:
% keys: a vector of structs that contain the keypoint and descriptors at
% all the various scales (blur sizes)
% 
% Author: Daniel Goodwin dgoodwin208@gmail.com 
% Date: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keys = SWITCH_tile_processingInParallel(img)
    
    loadExperimentParams;

    blur_size_list = params.SCALE_PYRAMID;

    keys_cell = cell(1,length(blur_size_list));

    parfor i = 1:length(blur_size_list)
        blur_size = blur_size_list(i);

        %Blurring is done inside the Harris keypoint detection code        
        res_vect = Harris3D(img, blur_size);

        %Blurring is done outside the 3D Sift code
        h  = fspecial3('gaussian',blur_size); 
        img_blur = convn(img,h,'same');         
        
        keys_cell{i} = calculate_3DSIFT(img_blur, res_vect);

    end
    
    keys = {};
    ctr = 1;

    for i = 1:length(blur_size_list)
        blur_size = blur_size_list(i);
        keys_blur = keys_cell{i};
        %Incrementally add each key idx, plus scale param k
        for key_idx=1:length(keys_blur)
            %sometimes the 3DSift code returns a 0 at the end. isnumeric checks for that
            if ~isnumeric(keys_blur{key_idx})
                %Note the level of the scale pyramid ()
                keys_blur{key_idx}.k = blur_size/2/2.354;
                keys{ctr} = keys_blur{key_idx};
                ctr = ctr+1;
            end
        end
    end

end
