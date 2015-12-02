function [ rgb ] = makeRGBImage(channel1, channel2, channel3 )
%MAKERGBIMAGE Summary of this function goes here
%   Detailed explanation goes here
% Now visualizing output image:


rgb = uint8(zeros(size(channel1,1),size(channel1,2),3));

r = 255*channel1/max(max(channel1));
rgb(:,:,1) = uint8(r);

if ~isempty(channel2) 
    g = 255*channel2/max(max(channel2))*.2;
    rgb(:,:,2) = uint8(g);
end

if ~isempty(channel3) 
    b = 255*channel3/max(max(channel3));
    rgb(:,:,3) = uint8(b);
end
end

