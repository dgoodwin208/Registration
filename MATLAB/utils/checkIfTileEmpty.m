function tileIsEmpty = checkIfTileEmpty(img,empty_threshold)

    total_xy_pix = size(img,1)*size(img,2);
    
    pix_averages = mean(img,3);
    
    num_black_pixels = sum(sum(pix_averages<1));
    
    percentage_empty = (num_black_pixels/total_xy_pix);
    if percentage_empty > empty_threshold
        tileIsEmpty = 1;
    else
        tileIsEmpty = 0;
    end
end
