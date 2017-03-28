function out = nearestInterpInParallel(imgw, map, maxhw )
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

[yi_arr,xi_arr,zi_arr] = ind2sub(size(map),I);

num_holes = length(yi_arr);

p = gcp('nocreate');
if isempty(p)
    stepsize = 1;
else
    stepsize = p.NumWorkers;
end

if isempty(yi_arr) == false

    yi_range = ceil(num_holes / stepsize);

    imgw_cell = cell(1,stepsize);
    map_cell  = cell(1,stepsize);

    % make slices of imgw and map
    for i = 1:stepsize
        ix_start = (i-1)*yi_range+1;
        ix_end   = i*yi_range;
        if i == stepsize
            ix_end = length(yi_arr);
        end

        disp(['i=' num2str(i) ' (' num2str(ix_start) ',' num2str(ix_end) ')'])

        yixL=max(min(yi_arr(ix_start:ix_end))-max(1,maxhw),1);
        yixU=min(max(yi_arr(ix_start:ix_end))+max(1,maxhw),outH);
        xixL=max(min(xi_arr(ix_start:ix_end))-max(1,maxhw),1);
        xixU=min(max(xi_arr(ix_start:ix_end))+max(1,maxhw),outW);
        zixL=max(min(zi_arr(ix_start:ix_end))-max(1,maxhw),1);
        zixU=min(max(zi_arr(ix_start:ix_end))+max(1,maxhw),outD);

        imgw_cell{i} = imgw(yixL:yixU,xixL:xixU,zixL:zixU);
        map_cell{i}  = map(yixL:yixU,xixL:xixU,zixL:zixU);

        disp(['wind=' num2str(yixL) ',' num2str(yixU) ';' num2str(xixL) ',' num2str(xixU) ';' num2str(zixL) ',' num2str(zixU)])
    end

    out_cell = cell(1,stepsize);

    disp('start parfor loop in nearestInterp')
    tic;
    disp(['stepsize=' num2str(stepsize) ', yi_range=' num2str(yi_range)])

    parfor i = 1:stepsize
        ix_start = (i-1)*yi_range+1;
        ix_end   = i*yi_range;
        if i == stepsize
            ix_end = length(yi_arr);
        end
        outH_sub = size(imgw_cell{i},1)
        outW_sub = size(imgw_cell{i},2)
        outD_sub = size(imgw_cell{i},3)

        disp(['i=' num2str(i) ' (' num2str(ix_start) ',' num2str(ix_end) ')'])
        out_cell{i} = zeros(1,ix_end-ix_start+1);

        yixBase=max(min(yi_arr(ix_start:ix_end))-max(1,maxhw),1);
        xixBase=max(min(xi_arr(ix_start:ix_end))-max(1,maxhw),1);
        zixBase=max(min(zi_arr(ix_start:ix_end))-max(1,maxhw),1);

        for ix = ix_start:ix_end
            
            xi = xi_arr(ix)-xixBase+1;
            yi = yi_arr(ix)-yixBase+1;
            zi = zi_arr(ix)-zixBase+1;

            % Find min window which has non-hole neighbors
            nz = false;
            yixL=0;
            yixU=0;
            xixL=0;
            xixU=0;
            zixL=0;
            zixU=0;
            for h = 1:maxhw
                yixL=max(yi-h,1);
                yixU=min(yi+h,outH_sub);
                xixL=max(xi-h,1);
                xixU=min(xi+h,outW_sub);
                zixL=max(zi-h,1);
                zixU=min(zi+h,outD_sub);
                
                if isempty(find(map_cell{i}(yixL:yixU,xixL:xixU,zixL:zixU), 1)) == false
                    nz = true;
                    break;
                end
            end

            % use the mean of non-hole neighbors in the window for filling
            if nz == true
                win = imgw_cell{i}(yixL:yixU, xixL:xixU, zixL:zixU);
                out_cell{i}(ix-ix_start+1) = mean(win(find(map_cell{i}(yixL:yixU, xixL:xixU,zixL:zixU)~=0)));
                if xi>outW
                    disp('ajh');
                end
                
            end
        end

    end

    disp('finished parfor loop in nearestInterp')
    toc;

    % gather each result from workers
    for i = 1:stepsize
        ix_start = (i-1)*yi_range+1;
        ix_end   = i*yi_range;
        if i == stepsize
            ix_end = length(yi_arr);
        end

        disp(['i=' num2str(i) ' (' num2str(ix_start) ',' num2str(ix_end) ')'])
        for ix = ix_start:ix_end
            if ~isempty(out_cell{i}(ix-ix_start+1))
                xi = xi_arr(ix);
                yi = yi_arr(ix);
                zi = zi_arr(ix);

                out(yi,xi,zi) = out_cell{i}(ix-ix_start+1);
            end
        end

    end

end

return;

