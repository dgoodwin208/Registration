#!/bin/bash

matlab -nosplash -nodisplay -r "profile = parallel.importProfile('scripts/cluster_profile__local_96workers'); parallel.defaultClusterProfile(profile); exit"

