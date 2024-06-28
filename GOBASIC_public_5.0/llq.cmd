#!/bin/ksh
# @ initialdir = /user/home/folder/GOBASIC/dir/
# @ requirements = (Arch == "ArchNumber") && (OpSys == "OSNumber")
# @ class = cluster_name
# @ notify_user = user@email
# @ error = error_file
# @ queue
#
# list commands to be executed below
cd $PBS_O_WORKDIR
./run_Run.sh /usr/local/MATLAB/R2014a > Run.out
