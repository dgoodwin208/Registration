function [m] = rotmat(a,b,c)
%calculates the rotation matrix in order XYZ in a right handed system
%the coordinate system is rotated, not the points!

ca=cos(a);
sa=sin(a);
cb=cos(b);
sb=sin(b);
cc=cos(c);
sc=sin(c);
m1=[1 0 0;0 ca -sa;0 sa ca];%rotation only about x axis
m2=[cb 0 sb;0 1 0;-sb 0 cb];%rotation only about y axis
m3=[cc -sc 0;sc cc 0;0 0 1];%rotation only about z axis
% m1=[1 0 0;0 ca sa;0 -sa ca];
% m2=[cb 0 -sb;0 1 0;sb 0 cb];
% m3=[cc sc 0;-sc cc 0;0 0 1];
m=m1*m2*m3;