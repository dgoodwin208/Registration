function  out_img = loadTifSequence( directory )

    %LOAD3DTIF: Load 3D a tif into a y,x,z stack
    %fname has to end in .tif

    files = dir(fullfile(directory,'*.tif'));
    
    num_images = length(files);
    
    A = imread(fullfile(directory,files(1).name));
    out_img = zeros([size(A),num_images]);
    
    sorted_names = sort_nat({files.name});
    for k = 1:num_images
        
        A = imread(fullfile(directory,sorted_names{k}));
        out_img(:,:,k)=A;
    end

end

