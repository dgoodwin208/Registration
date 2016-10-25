function  out_img = loadTifSequence( directory,start_idx,end_idx )

    %LOAD3DTIF: Load 3D a tif into a y,x,z stack
    %fname has to end in .tif
    
    files = dir(fullfile(directory,'*.tif'));
    
    num_images = length(files);
    
    A = imread(fullfile(directory,files(1).name));
    
    sorted_names = sort_nat({files.name});
    
    %If the user didn't specify file start index, start at the beginning
    if ~exist('start_idx','var')
        start_idx = 1;
    end
    %Similar for default loading everything
    if ~exist('end_idx','var')
        end_idx = num_images;
    else
        %Don't let the user specify more files than there are available
        end_idx = min(num_images,end_idx);
    end
    
    
    out_img = zeros([size(A),end_idx-start_idx+1]);
    
    idx = 1;
    for k = start_idx:end_idx
        
        A = imread(fullfile(directory,sorted_names{k}));
        out_img(:,:,idx)=A;
        idx = idx+1;
    end

end

