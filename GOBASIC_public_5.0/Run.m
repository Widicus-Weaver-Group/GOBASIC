%clear
tic_id= tic;
definevars;
analysisscript;
anonfun4;
toc(tic_id)
% Enable at parallel computing cluster
% fitjob = parpool('local', 32); % input the number of cores on your cluster
optest;
% delete(fitjob)
%exit
