function saveTifSequence(imgvol, outputDirectory)

    %SAVE3DTif: Save a 3D volume into an image sequence. To be used in the
    %case of the vol being too big to be loaded as a single 3D tif
    
    
    %set the minimum to zero
    imgvol = imgvol - min(min(min(imgvol)));
    
    %set the max to max uint16
    imgvol = (imgvol/max(max(max(imgvol))))*(2^16-1);
    imgvol = uint16(imgvol);
    
    if exist(outputDirectory) ~= 7
        mkdir(outputDirectory);
    end
    
    
   for K=1:size(imgvol,3)
    outputFileName = fullfile(outputDirectory, [num2str(K),'.tif']);
    imwrite(imgvol(:, :, K), outputFileName, 'Compression','none');
   end
    
    
end
