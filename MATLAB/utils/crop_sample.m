%%Crop each channel based on the xtrema of non-emmpty pixels based on
%%Lectin channels of each sample
function crop_sample(sample_num, run_num)

    loadExperimentParams;
    
    filename = fullfile(params.INPUTDIR,sprintf('sample%dround%d_%s.tif',sample_num,run_num,'Lectin')); 
    disp('Loading Data...');
    img = load3DTif(filename);
    pix_averages = mean(img,3);

    disp('Generating stacked img based on pixel z-dim mean');
    disp('Starting MATLABs imcrop. Draw the bounding box then double click');
    [I2, rect] = imcrop(pix_averages);

    imagesc(img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),25))
    
    y_min = round(rect(2));
    x_min = round(rect(1));
    y_max = round(rect(2) + rect(4));
    x_max = round(rect(1) + rect(3));
    
    bounds = [y_min, y_max, x_min, x_max];
    cropfilename = fullfile(params.OUTPUTDIR,sprintf('sample%dround%d_cropbounds.mat',params.SAMPLE_NUM,run_num));
    save(cropfilename,'bounds')
    disp(['Crop bounds saved to ' cropfilename]);
end