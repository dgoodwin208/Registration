%If you have the MATLAB parallel computing toolbox, a script like this 
%may be useful for running the time intensive, yet embarassingly parallel
%process of calculating keypoints and descriptors in each 

%without specifying, MATLAB uses 
parpool(); 
loadExperimentParams;
round_number_to_process = 2;

parfor i = 7:params.ROWS_DESC*params.COLS_DESC
    calculateDescriptors(round_number_to_process,i,i);
end

poolobj = gcp('nocreate');
delete(poolobj);
