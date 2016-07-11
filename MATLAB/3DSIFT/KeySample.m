function index = KeySample(key, pix)

LoadParams;

fv = sphere_tri('ico',Tessellation_levels,1);


irow = int16(key.x);
icol = int16(key.y);
islice = int16(key.z);

xySpacing = key.xyScale * MagFactor;
tSpacing = key.tScale * MagFactor;

xyRadius = 1.414 * xySpacing * (IndexSize + 1) / 2.0;
tRadius = 1.414 * tSpacing * (IndexSize + 1) / 2.0;
xyiradius = int16(xyRadius);
tiradius = int16(tRadius);

index = zeros(IndexSize,IndexSize,IndexSize,nFaces);

for i = -xyiradius:xyiradius
    for j = -xyiradius:xyiradius
        for s = -tiradius:tiradius

            % This is redundant and probably slows down the code, but at
            % some point this solved a major overflow headache, so leaving
            % as-is for now
            i2 = double(i);
            j2 = double(j);
            s2 = double(s);
            distsq = double(i2^2 + j2^2 + s2^2);

            v0 = [i2; j2; s2];

            
            i_indx = int16(floor(double((i + xyiradius)) / double((2*xyiradius/IndexSize)))) + 1;
            j_indx = int16(floor(double((j + xyiradius)) / double((2*xyiradius/IndexSize)))) + 1;
            s_indx = int16(floor(double((s + tiradius)) / double((2*tiradius/IndexSize)))) + 1;
            
            if i_indx > IndexSize
                i_indx = IndexSize;
            end
            if j_indx > IndexSize
                j_indx = IndexSize;
            end
            if s_indx > IndexSize
                s_indx = IndexSize;
            end

            if (i_indx < 1 || j_indx < 1 || s_indx < 1)
                disp('Something wrong with the sub-histogram index');
            end
            
            %For each pixel, take a neighborhhod of xyradius and tiradius,
            %bin it down to the IndexSize dimensions
            r = irow + v0(1);
            c = icol + v0(2);
            t = islice + v0(3);

            %We've calculated the target index for the 3D histogram
            %Add sample increments the 
            index = AddSample(index, pix, distsq, r, c, t, i_indx, j_indx, s_indx, fv);
            
        end
    end
end

%IF NOT ADDING ROT INVARIANCE, RETURN HERE.
return;

% disp('DOING ROT INV STUFF');
%For each tesselation face, add up the contribution from each face
sum_mags = size(index,4);
for t = 1:size(index,4)
    sum_mags(t) = sum(sum(sum(index(:,:,:,t))));
end



[mags_sorted, mags_indexed] = sort(sum_mags,'descend');

%Get the most dominant face:
vect_dominant = fv.centers(mags_indexed(1),:);

theta = atan(vect_dominant(1)/vect_dominant(2));
phi = atan(vect_dominant(3)/sqrt(vect_dominant(1)*vect_dominant(1) + vect_dominant(2)*vect_dominant(2)));

rotate_matrix = [cos(theta)*cos(phi), -1*sin(theta),-1*cos(theta)*sin(phi), 0;
                 sin(theta)*cos(phi), cos(theta),   -1*sin(theta)*sin(phi), 0;
                 sin(phi),            0,            cos(phi),               0;
                 0,                   0,            0,                      1];

             
             
%--------------------------------------------------             
%Rotating tesseltated sphere approach             
rotated = [fv.centers, ones(size(fv.centers,1),1)]*rotate_matrix';
rotated = rotated(:,1:3); %There was really no good reason to keep the homogenous coords

%For each rotated point, find the index of the original face
%this will be then used to reorient the histogram
re_indexed = zeros(size(fv.centers,1),1);

%We want to find the index to map from the original space (fv.centers) to
%the new rotated space
for j = 1:size(fv.centers,1)
    %initialize the pointers for the finding the nearest original position
    min_dist = 1000;
    min_idx = -1;
    
    %loop over all original locations to see which the new one is closest
    for i=1:size(rotated,1)
        dist = sqrt(  (rotated(i,1) - fv.centers(j,1))^2 + (rotated(i,2) - fv.centers(j,2))^2 + (rotated(i,3) - fv.centers(j,3))^2 );
        if dist<min_dist
%             disp([mat2str(rotated(i,:)) ' ' mat2str(fv.centers(j,:)) ' ' num2str(dist)]);
            min_dist = dist;
           min_idx = i;
        end
    end
    %re_indexed is the link from the original space to rotated
    re_indexed(j) = min_idx;
end

% DEBUGGING TOOL
% new_mapped_to_old = fv.centers(re_indexed,:);
% 
% figure
% scatter3(fv.centers(:,1), fv.centers(:,2), fv.centers(:,3))
% hold on;
% 
% %scatter3(rotated(:,1), rotated(:,2), rotated(:,3),'r')
% for j = 1:size(fv.centers,1)
%    point_pair = [fv.centers(j,:);new_mapped_to_old(j,:)];
%    scatter3(new_mapped_to_old(j,1), new_mapped_to_old(j,2), new_mapped_to_old(j,3),'r')
%    line(point_pair(:,1), point_pair(:,2), point_pair(:,3),'Color',[0 0 0]);
% end
% hold off;


%Moving the original histogram space into the rotated space requires one
%last sum
rotated_index = zeros(IndexSize,IndexSize,IndexSize,nFaces);
for j = 1:size(fv.centers)
    rotated_index(:,:,:,j) = rotated_index(:,:,:,j) + index(:,:,:,re_indexed(j));
end

% figure;
% plot(index(:));hold on;
% plot(rotated_index(:)); hold off;

index = rotated_index; 

return;             
             
             
             
             
             
             
             
             
%--------------------------------------------------             
%Rotating neighborhood approach:
%Now, we can either rotate the whole image and recalcuate the features (stored in index)
%or, we can recalculatd index around a rotated_pix object. 
%As we hit issues trying to rotate the tesselated sphere, let's try
%rotating pix around (irow, icol, islice)

% t = [irow; icol; islice] - rotate_matrix(1:3,1:3)*[irow; icol; islice];
% rotate_matrix(1:3,4) = t;

%The biggest rotation that would roatate the matrix is 45, which would
%increase the size of the image by sqrt(2)
padfactorxy = ceil(sqrt(2)*xyRadius);
padfactort = ceil(sqrt(2)*tiradius);

%Make a subregion around the keypoint
subpix = pix(irow-padfactorxy:irow+padfactorxy, icol-padfactorxy:icol+padfactorxy, islice - padfactort:islice+padfactort);

%Create a transform matrix to move to the center pixel (which we'll use
%twice)

%recalculate the position of the keypoint in the new mini space
irow = int16(padfactorxy+1);
icol = int16(padfactorxy+1);
islice = int16(padfactort+1);

translate_matrix_fwd = double([eye(4,3), [-irow; -icol; -islice; 1]]);
translate_matrix_back = double([eye(4,3), [irow; icol; islice; 1]]);

%Rotate that subregion
tformObj = affine3d(translate_matrix_fwd' * rotate_matrix' * translate_matrix_back');
Rcb = imref3d(size(subpix));
[outputImage,rb] = imwarp(subpix,tformObj,'OutputView',Rcb, 'Interp', 'cubic');

% figure(1);
% subplot(1,2,1)
% imagesc(subpix(:,:,islice));
% subplot(1,2,2)
% subpix = outputImage;
% imagesc(subpix(:,:,islice));



%Now re-run the calculation of the index
index = zeros(IndexSize,IndexSize,IndexSize,nFaces);
for i = -xyiradius:xyiradius
    for j = -xyiradius:xyiradius
        for s = -tiradius:tiradius

            % This is redundant and probably slows down the code, but at
            % some point this solved a major overflow headache, so leaving
            % as-is for now
            i2 = double(i);
            j2 = double(j);
            s2 = double(s);
            distsq = double(i2^2 + j2^2 + s2^2);

            v0 = [i2; j2; s2];

            
            i_indx = int16(floor(double((i + xyiradius)) / double((2*xyiradius/IndexSize)))) + 1;
            j_indx = int16(floor(double((j + xyiradius)) / double((2*xyiradius/IndexSize)))) + 1;
            s_indx = int16(floor(double((s + tiradius)) / double((2*tiradius/IndexSize)))) + 1;
            
            if i_indx > IndexSize
                i_indx = IndexSize;
            end
            if j_indx > IndexSize
                j_indx = IndexSize;
            end
            if s_indx > IndexSize
                s_indx = IndexSize;
            end

            if (i_indx < 1 || j_indx < 1 || s_indx < 1)
                disp('Something wrong with the sub-histogram index');
            end
            
            %For each pixel, take a neighborhhod of xyradius and tiradius,
            %bin it down to the IndexSize dimensions
            r = irow + v0(1);
            c = icol + v0(2);
            t = islice + v0(3);

            %We've calculated the target index for the 3D histogram
            %Add sample increments the 
            %Even though key is a parameter (which is now in a differnt
            %coord space, it's not used in the AddSample function)
            index = AddSample(index, subpix, distsq, r, c, t, i_indx, j_indx, s_indx, fv);
            
        end
    end
end

return
             
%              
% rotated = [fv.centers, ones(size(fv.centers,1),1)]*rotate_matrix';
% rotated = rotated(:,1:3); %There was really no good reason to keep the homogenous coords
% 
% %For each rotated point, find the index of the original face
% %this will be then used to reorient the histogram
% re_indexed = zeros(size(fv.centers,1),1);
% 
% %We want to find the index to map from the original space (fv.centers) to
% %the new rotated space
% for j = 1:size(fv.centers,1)
%     %initialize the pointers for the finding the nearest original position
%     min_dist = 1000;
%     min_idx = -1;
%     
%     %loop over all original locations to see which the new one is closest
%     for i=1:size(rotated,1)
%         dist = sqrt(  (rotated(i,1) - fv.centers(j,1))^2 + (rotated(i,2) - fv.centers(j,2))^2 + (rotated(i,3) - fv.centers(j,3))^2 );
%         if dist<min_dist
% %             disp([mat2str(rotated(i,:)) ' ' mat2str(fv.centers(j,:)) ' ' num2str(dist)]);
%             min_dist = dist;
%            min_idx = i;
%         end
%     end
%     %re_indexed is the link from the original space to rotated
%     re_indexed(j) = min_idx;
% end
% 
% % DEBUGGING TOOL
% % new_mapped_to_old = fv.centers(re_indexed,:);
% % 
% % figure
% % scatter3(fv.centers(:,1), fv.centers(:,2), fv.centers(:,3))
% % hold on;
% % 
% % %scatter3(rotated(:,1), rotated(:,2), rotated(:,3),'r')
% % for j = 1:size(fv.centers,1)
% %    point_pair = [fv.centers(j,:);new_mapped_to_old(j,:)];
% %    scatter3(new_mapped_to_old(j,1), new_mapped_to_old(j,2), new_mapped_to_old(j,3),'r')
% %    line(point_pair(:,1), point_pair(:,2), point_pair(:,3),'Color',[0 0 0]);
% % end
% % hold off;
% 
% %index = index(:,:,:,re_indexed);
% %Moving the original histogram space into the rotated space requires one
% %last sum
% rotated_index = zeros(IndexSize,IndexSize,IndexSize,nFaces);
% for j = 1:size(fv.centers)
%     rotated_index(:,:,:,j) = rotated_index(:,:,:,j) + index(:,:,:,re_indexed(j));
% end
% 
% index = rotated_index;
