set -e #exit the script if there are errors; this guarantees that both syncs were successful if the script runs through

rclone copy /Users/dp/projects/RVA_2 Duke_box:/PROJECT\ 3373\:\ Projects/RVA_2 -v #copies to Box
