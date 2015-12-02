%If you have the MATLAB parallel computing toolbox, a script like this 
%may be useful for running the time intensive, yet embarassingly parallel
%process of calculating keypoints and descriptors in each 

%without specifying, MATLAB uses 
parpool(); 
loadExperimentParams;
experiment_number_to_process = 1;
parfor i = 1:params.ROWS_DESC*params.COLS_DESC
    CalculateDescriptorsForTileAtIndices(experiment_number_to_process ,i,i);
end