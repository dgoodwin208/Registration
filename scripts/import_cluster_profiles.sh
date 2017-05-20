#!/bin/bash

matlab -nosplash -nodisplay -r "p=parallel.clusterProfiles; if ~find(ismember(p,'local_96workers')); profile = parallel.importProfile('tests/cluster_profile__local_96workers'); parallel.defaultClusterProfile(profile); else disp('already imported.'); end; exit"

