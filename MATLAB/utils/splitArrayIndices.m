function [indices] = splitArrayIndices(inputarray_dims, dim, chunksize)
% Similar to linspace, given an array and a desired chunk size, this
% function returns an array of indices separated by chunksize. The last
% entry may be less than chunksize.
indices = [1]; ctr = 1;

while (ctr<inputarray_dims(dim))
    next = min(indices(end)+chunksize,inputarray_dims(dim));
    
    
    %If the gap between the last two indices is less than half the desired
    %chunk size, then just merge the last chunk in with 
    %consider this edge case: (just merge the last two sections
    % 1    51   101   151   201   251   301   351   401   451   453
    if next - indices(end) < .5*chunksize
        indices(end) = next;
    else
        indices(end+1) = next;
    end
    ctr = ctr+chunksize;
end

end
