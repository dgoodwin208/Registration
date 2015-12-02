function [ outputImage_interp ] = TPS3DApply(in1D,out1D,imgToWarp,outdim)
%TPS3DAPPLY applies the warp to a vectorized version of the image data,
%then returns an interpolated version of the data

                          
%now that we have the pixel->pixel mapping from the TPS, apply it to a 
%data channel, interpolate, remove the padding, then save to disk
nominal_offset = .0001; %used to find holes of TPS Warp
outputImage = zeros(outdim);


%do the mapping from the input to the output image
%add the nominal offset to make all mapped pixels non-zero
try 
    outputImage(out1D) = imgToWarp(in1D)+nominal_offset;
catch
    disp('WARNING: There is an indexing error in the TPS Application ');    
end

%interpolate
outputImage_interp = interpolateVolume(outputImage);

%remove the nominal offset
outputImage_interp = round(outputImage_interp); 

end

