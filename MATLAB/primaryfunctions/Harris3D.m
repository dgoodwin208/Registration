% Step 4: Keypoints and Descriptors
% This is the code that calculates the keypoints locations via Harris Corner 3D
%
% INPUTS: 
% img is the image volume
% blur_dim describes the size of gaussian filter
% 
% OUTPUTS:
% res_vect: an n x 3 vector of coordinates of distinctive keypoints 
% 
% Author: Daniel Goodwin dgoodwin208@gmail.com 
% Date: August 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  res_vect = Harris3D(test_img,blur_dim)

%Creating 3D sobel filters
h = [ -1; 2; 1];
hp = [1 0 -1];

Hx = zeros(3,3,3);
Hy = zeros(3,3,3);
Hz = zeros(3,3,3);

Hy(1,:,:) = hp(1)*h*h';
Hy(2,:,:) = hp(2)*h*h';
Hy(3,:,:) = hp(3)*h*h';

Hx(:,1,:) = hp(1)*h*h';
Hx(:,2,:) = hp(2)*h*h';
Hx(:,3,:) = hp(3)*h*h';

Hz(:,:,1) = hp(1)*h*h';
Hz(:,:,2) = hp(2)*h*h';
Hz(:,:,3) = hp(3)*h*h';



% making filters (note that we're ignoring the fx, fy fz now)
clear fx fy fz

% applying sobel edge detector in the horizontal direction
Ix = convn(test_img,Hx,'same');
disp('calc sobel of x dir')

Iy = convn(test_img,Hy,'same'); 
disp('calc sobel of y dir')

Iz = convn(test_img,Hz,'same');
disp('calc sobel of z dir')

Ix2 = Ix.^2;
Iy2 = Iy.^2;
Iz2 = Iz.^2;
Ixy = Ix.*Iy;
Ixz = Ix.*Iz;
Iyz = Iy.*Iz;



tic
disp('calculating gaussian filters')

% applying gaussian filter on the computed value
h   = fspecial3('gaussian',blur_dim);

Ix2 = convn(Ix2,h,'same'); 
Iy2 = convn(Iy2,h,'same');
Iz2 = convn(Iz2,h,'same');

Ixy = convn(Ixy,h,'same');
Ixz = convn(Ixz,h,'same');
Iyz = convn(Iyz,h,'same');

height = size(test_img,1);
width = size(test_img,2);
depth = size(test_img,3);
result = zeros(height,width,depth); 
%R is the 'Corner response' matrix as defined in the original 1988 paper
R = zeros(height,width,depth);
toc 

tic
disp('running getting Rmax')


kappa = 1/27; %Magic number from the paper
margin = floor(.1*size(test_img,3));

% det(M) - kappa*trace(M)^2
R = Ix2.*(Iy2.*Iz2-Iyz.*Iyz)-Ixy.*(Ixy.*Iz2-Ixz.*Iyz)+Ixz.*(Ixy.*Iyz-Ixz.*Iy2)-kappa*(Ix2+Iy2+Iz2).^2;

% Clip R with margin
R(1:margin-1,:,:) = 0;
R(height-margin+1:height,:,:) = 0;
R(:,1:margin-1,:) = 0;
R(:,width-margin+1:width,:) = 0;
R(:,:,1:margin-1) = 0;
R(:,:,depth-margin+1:depth) = 0;


vecsort = sort(R(:),'descend');

%Because of the strong power law in this data, choose RMax to be the value of the 50th highest
Rmax = vecsort(5000); 

toc

%% Generating the points

 %If the candidate spot is of high enough value
[r,c,v] = ind2sub(size(R),find(R > .1*Rmax));
cnt = 0;
disp('Looping over candidate positions')
res_vect=[];

%mm = 5;% max margin: the neighborhood around which we get a maximum
mm = ceil((blur_dim/2/2.354));
%For a cube of 2*mm+1 per side, what is the index of the middle voxel?

center_idx = sub2ind([2*mm+1,2*mm+1,2*mm+1],mm+1,mm+1,mm+1);

for ctr=1:length(r)
    i = r(ctr);
    j = c(ctr);
    k = v(ctr);
    
    %avid boundary problems
    if i<=mm || j <=mm || k <=mm || i>=height-mm+1 || j>=width-mm+1 || k>=depth-mm+1
        continue
    end
    
    nbd = R(i-mm:i+mm,j-mm:j+mm,k-mm:k+mm);
    idx = find(nbd==max((max(max(nbd)))));

    %Take it as the max of the region if it's in the center
    if length(idx)==1
        if idx==center_idx
            result(i,j,k) = 1;
            cnt = cnt+1;
            res_vect(:,cnt) = [i j k];
        end
    else %If there is a tie
        disp('Seeing a tie')
        for index=1:length(idx)
            if idx(index)==center_idx
                result(i,j,k) = 1;
                res_vect(:,cnt) = [i j k];
                cnt = cnt+1;
            end
        end

    end;

end;

res_vect = res_vect';

end




