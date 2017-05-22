#!/bin/bash

matlab -nosplash -nodisplay -r "p=parallel.clusterProfiles; if sum(ismember(p,'local_96workers'))==0; profile = parallel.importProfile('scripts/cluster_profile__local_96workers'); parallel.defaultClusterProfile(profile); else disp('already imported.'); end; exit"

