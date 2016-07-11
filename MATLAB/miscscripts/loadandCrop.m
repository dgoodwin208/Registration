
disp(['RUNNING ON MOVING: ' num2str(params.MOVING_RUN) ', FIXED: ' num2str(params.FIXED_RUN) ', ' params.DATACHANNEL])

filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
                    params.SAMPLE_NAME,params.FIXED_RUN,params.DATACHANNEL));

imgFixed_total = load3DTif(filename);


filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
                    params.SAMPLE_NAME,params.MOVING_RUN,params.DATACHANNEL));

imgMoving_total = load3DTif(filename);

%LOAD FILES WITH CROP INFORMATION, CROP LOADED FILES
cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%d_cropbounds.mat',params.SAMPLE_NAME,params.FIXED_RUN));
if exist(cropfilename,'file')==2
    load(cropfilename,'bounds'); bounds_fixed = floor(bounds); clear bounds;
    imgFixed_total = imgFixed_total(bounds_fixed(1):bounds_fixed(2),bounds_fixed(3):bounds_fixed(4),:);        
end

cropfilename = fullfile(params.OUTPUTDIR,sprintf('%sround%d_cropbounds.mat',params.SAMPLE_NAME,params.MOVING_RUN));
if exist(cropfilename,'file')==2
    load(cropfilename,'bounds'); bounds_moving = floor(bounds); clear bounds;
    imgMoving_total = imgMoving_total(bounds_moving(1):bounds_moving(2),bounds_moving(3):bounds_moving(4),:);
end