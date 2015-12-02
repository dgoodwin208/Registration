function  out_img = load3DTif( fname )

    %LOAD3DTIF: Load 3D a tif into a y,x,z stack
    %fname has to end in .tif

    info = imfinfo(fname);
    num_images = numel(info);
    A = imread(fname, 1, 'Info', info);
    out_img = zeros([size(A),num_images]);
    for k = 1:num_images
        A = imread(fname, k, 'Info', info);
        out_img(:,:,k)=A;
    end

end

