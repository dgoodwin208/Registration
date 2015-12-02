function out = nearestInterp(imgw, map, maxhw )
% Description:
% Fill holes using nearest neighbor (median filter of neighbors) interpolation
%
% Inputs:
% imgw - input image
% map - Map of the canvas with 0 indicating holes and 1 indicating pixel
% maxhw - Max radius for nearest neighbor interpolation
%
% Output:
% out - interpolated image
%
% Author: Fitzgerald J Archibald
% Date: 23-Apr-09

outH  = size(imgw,1);
outW  = size(imgw,2);
outD  = size(imgw,3);
out = imgw;

I = find(map==0); % Find locations needing fill

[yi_arr, xi_arr,zi_arr] = ind2sub(size(map),I);

num_holes = length(yi_arr);

ctr=1;
tic; %start the timer.
if isempty(yi_arr) == false
    
    for ix = 1:length(yi_arr),
        
        xi = xi_arr(ix);
        yi = yi_arr(ix);
        zi = zi_arr(ix);

        % Find min window which has non-hole neighbors
        nz = false;
        for h = 1:maxhw,
            yixL=max(yi-h,1);
            yixU=min(yi+h,outH);
            xixL=max(xi-h,1);
            xixU=min(xi+h,outW);
            zixL=max(zi-h,1);
            zixU=min(zi+h,outD);
            
            if isempty(find(map(yixL:yixU,xixL:xixU,zixL:zixU), 1)) == false
                nz = true;
                break;
            end
        end

        % use the mean of non-hole neighbors in the window for filling
        if nz == true,
            win = imgw(yixL:yixU, xixL:xixU, zixL:zixU);
            out(yi,xi,zi) = mean(win(find(map(yixL:yixU, xixL:xixU,zixL:zixU)~=0)));
            if xi>outW
                disp('ajh');
            end
            
        end
        
        ctr= ctr+1;
        %run a time estimate once
        
        if mod(ctr,100000)==0 
           time_elapsed = toc;
           disp(['Estimated remaining time for interpolation: ' num2str(time_elapsed*(num_holes-ctr)/100000)]);
           tic; %restart the timer
        end
    end
end

return;

