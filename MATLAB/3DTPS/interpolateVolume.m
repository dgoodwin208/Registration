function outputImage = interpolateVolume(inputImage)
%Assumes that there are 0s in the input image where the original image was
%not able to warp into. To ensure this, a nominal f

%Add a nominal amount to the tileM for the interpolation step at the end of
%this function
%It is rounded back to integers at the very end.


%make a max_projection 2D image then make a mask of it
se = strel('disk',10); %for the morphological closing operation
mask2d = imclose(max(inputImage,[],3),se);
mask_indices = find(mask2d>0);
[mask_y,mask_x] = ind2sub([size(inputImage,1),size(inputImage,2)], mask_indices);

%Here we create a map of the pixels that were not set due to the rounding
%errors of the TPS.
%Method is we zoom into a tight rectangular crop around the data, then use
%the imclosed version of the maxintensity projection to search for any
%pixels that are zero inside the mask (by adding the nominal amount only
%unset pixels should be zero)
rect_section = inputImage(min(mask_y):max(mask_y), min(mask_x):max(mask_x), :);
map = ones(size(rect_section));

mask2d_small = imclose(max(rect_section,[],3),se);
mask2d_small_mask = mask2d_small>0.1; %what indices are in the closed version of max projection
for z_idx = 1:size(rect_section,3)
    
    %Creating a binary volume to tell the nearestInterp function where the
    %values are that need interpolation (0 indicates a hole)
    slice_out = ones([size(rect_section,1),size(rect_section,2)]);
    slice_data = rect_section(:,:,z_idx);
    slice_out(mask2d_small_mask) = slice_data(mask2d_small_mask)>0;
    map(:,:,z_idx) = slice_out;    

end

outputImage = inputImage;

%Empirically, interp_radius of 1 creates some artifacts (mainly visible
%around areas of high contrast). interp_radius of 2 is slower (within one order of magnitude) 
%but has better results
interp_radius = 2;
interpolated_rect = nearestInterpInParallel(rect_section, map, interp_radius );
%interpolated_rect = nearestInterp(rect_section, map, interp_radius );

outputImage(min(mask_y):max(mask_y), min(mask_x):max(mask_x), :) = interpolated_rect ;

end
