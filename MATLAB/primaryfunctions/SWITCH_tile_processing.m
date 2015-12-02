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

function keys = SWITCH_tile_processing(img)
    
    keys_total = {};
    ctr = 1;

    loadExperimentParams;
    
    for blur_size= params.SCALE_PYRAMID
        
        %Blurring is done inside the Harris keypoint detection code        
        res_vect = Harris3D(img, blur_size);

        %Blurring is done outside the 3D Sift code
        h  = fspecial3('gaussian',blur_size); 
        img_blur = convn(img,h,'same');         
        
        keys = calculate_3DSIFT(img_blur, res_vect);
        
        %Incrementally add each key idx, plus scale param k
        for key_idx=1:length(keys)
            %sometimes the 3DSift code returns a 0 at the end. isnumeric checks for that
            if ~isnumeric(keys{key_idx})
                %Note the level of the scale pyramid ()
                keys{key_idx}.k = blur_size/2/2.354;
                keys_total{ctr} = keys{key_idx};
                ctr = ctr+1;
            end
        end
    end
    
    keys = keys_total;

end