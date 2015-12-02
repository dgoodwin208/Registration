% Code taken from Paul Scovanner's homepage: 
% http://www.cs.ucf.edu/~pscovann/


% Calculate 3D Sift
% img is the input in question.
% keypts is an Nx3 matrix of keypoints

function  keys = calculate_3DSIFT(img, keypts)

i = 0;
offset = 0;
while 1

    reRun = 1;
    i = i+1;
    
    while reRun == 1
        
        loc = keypts(i+offset,:);
        %fprintf(1,'Calculating keypoint at location (%d, %d, %d)\n',loc);
        
        % Create a 3DSIFT descriptor at the given location
        [keys{i} reRun] = Create_Descriptor(img,1,1,loc(1),loc(2),loc(3));
        
        if reRun == 1
            offset = offset + 1;
        end
        
        %are we out of data?
        if i+offset>=size(keypts,1)
            break;
        end
    end
    
    %are we out of data?
    if i+offset>=size(keypts,1)
            break;
    end
end

fprintf(1,'\nFinished.\n%d points thrown out do to poor descriptive ability.\n',offset);

end