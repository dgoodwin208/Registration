function [ gridimg ] = makeGrid( dims )
%MAKEGRID Make a white cube grid of lines to show the effects of the TPS
%   dims are the shape of the 3D tif that you want to mimic

gridimg =  zeros(dims);
xyspacing = 100;
zspacing = 10;
%get the vertical lines
indices = 1:size(gridimg,1);
vI = find(mod(indices,xyspacing)==1);
for z=1:zspacing:size(gridimg,3)
    gridimg(vI,:,z:z+2) = 255;
    gridimg(vI+1,:,z:z+2) = 255;
    gridimg(vI+2,:,z:z+2) = 255;
end

%get the horizontal lines
indices = 1:size(gridimg,2);
hI = find(mod(indices,xyspacing)==1);

for z=1:zspacing:size(gridimg,3)
    gridimg(:,hI,z:z+2) = 255;
    gridimg(:,hI+1,z:z+2) = 255;
    gridimg(:,hI+2,z:z+2) = 255;    
end

%draw the depth lines at the intersection of vertical and horizontal lines
for y=vI
    for x=hI
        gridimg(y:y+2,x:x+2,:) = 255;
%         gridimg(y+1,x,:) = 255;
%         gridimg(y+1,x+1,:) = 255;                
%         gridimg(y,x+1,:) = 255;                
    end
end

    
% [VI,HI] = meshgrid(vI,hI);
% gridimg(VI,HI,:) = 255;
% gridimg(VI+1,HI+1,:) = 255;

end

