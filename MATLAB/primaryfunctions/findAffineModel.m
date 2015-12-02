function  affine_tform = findAffineModel(key1, key2)

num_samples = size(key1,1);

%Put the coordinates into homogenous coords (appending the 1)   
x1 = [key1, ones(num_samples,1)];

%Format the data into a matrix that we will take the pseudo inverse of
A = zeros(num_samples*3, 12);
for i = 1:num_samples
    A(3*(i-1)+1:3*i, :) = [ x1(i,:) zeros(1,8); zeros(1,4), x1(i,:), zeros(1,4); zeros(1,8), x1(i,:)]; 
end


%Vectorize the output points (x1 mapped into x2)
x2 = key2';
b = x2(:);

%debug sanity check printouts
% key1(1:4,:)
% key2(1:4,:)
% A(1:12,:)
% b(1:13)

%Calculate the params to go from space 1 to space 2
params = A\b;
%Condition matrix can be checked with cond(A);
affine_tform = [params(1) params(2) params(3) params(4);
                params(5) params(6) params(7) params(8);
                params(9) params(10) params(11) params(12);
                0           0       0           1; ];


            
end