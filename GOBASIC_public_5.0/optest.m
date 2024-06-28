% *************************************************************************
% Description: Optimization and output
% Author: Mary L. Radhuber
% Date: May 2, 2013
% Developer: Luyao Zou
% v2.0: Dec. 25, 2013
% v3.0: May. 29, 2014
% v4.0: Sept. 13, 2014
% v4.3: Feb. 15, 2014
% v4.4: Jun. 17, 2015
% *************************************************************************

% >>>>>>>>>>>> FIT SECTION >>>>>>>>>>>> 
tic_id = tic; 
[x,fval] = patternsearch(fitfunc,init,[],[],[],[],low,up,[],options);
toc(tic_id);

% >>>>>>>>>>>> OUTPUT PREPERATION SECTION >>>>>>>>>>>> 
% 
fits = spect(x);
totalfit = sum(fits,2) + update;
residual = (totalfit - obsint).*nonflagidx;

% Determine ta_err from observational data OR user input
if ta_err_from_fit
    ta_err = norm(residual)/sqrt(length(residual)-size(x,1)*comp_num);
else
    ta_err = ta_err_deconv;
end

% Calculate uncertainty of parameters
% uncertainty rows: [u_log10NT, u_dv, u_0.01temp, u_shift, u_NT]'
tic_id = tic;
fprintf('\nEstimating Uncertainties\n');
uncertainty = err(x,tol,residual,ta_err,molfreq,obsfreq,aco,agu,...
                  eup,lc_idx,chwidth,beamfilling);
toc(tic_id);

% Convert log10NT & 0.01temp back in row 1 & 3
% results rows: [log10NT, dv, Temp, shift, NT]'
results = [x(1,:); uncertainty(1,:);...
           x(2,:); uncertainty(2,:);...
           100.*x(3,:); 100.*uncertainty(3,:);...
           x(4,:); uncertainty(4,:);...
           10.^x(1,:); 10.^x(1,:).*log(10).*uncertainty(1,:);];

% >>>>>>>>>>>> OUTPUT SECTION >>>>>>>>>>>> 
% Display results on the Matlab terminal and make plots
displayresult;
% Write results into files
writeresult;
