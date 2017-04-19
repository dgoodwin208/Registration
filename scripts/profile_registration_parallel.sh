#!/bin/bash

### calculateDescriptors for round 1 and 2 in parallel
matlab -nodisplay -nosplash -logfile ./matlab-calcDesc.log -r "addpath(genpath('./MATLAB')); calculateDescriptorsInParallel([1 2]); exit"

#### registerWithDescriptors for round 1 and 2
matlab -nodisplay -nosplash -logfile ./matlab-registerWDesc.log -r "addpath(genpath('./MATLAB')); registerWithDescriptors(2); exit"

exit 0

