#!/bin/bash

### calculateDescriptors for round 1
matlab -nodisplay -nosplash -logfile ./matlab-calcDesc1.log -r "addpath(genpath('./MATLAB')); calculateDescriptors(1); exit"

### calculateDescriptors for round 2
matlab -nodisplay -nosplash -logfile ./matlab-calcDesc2.log -r "addpath(genpath('./MATLAB')); calculateDescriptors(2); exit"

### registerWithDescriptors for round 1 and 2
matlab -nodisplay -nosplash -logfile ./matlab-registerWDesc.log -r "addpath(genpath('./MATLAB')); registerWithDescriptors(2); exit"

exit 0

