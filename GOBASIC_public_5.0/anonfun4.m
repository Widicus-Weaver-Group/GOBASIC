% *************************************************************************
% Description: Include previous fits and generate fitting function.
% Author: Mary Radhuber
% Date : May 2, 2013
% Developer: Luyao Zou
% v2.0: Dec. 25, 2013
% v2.0: Jan. 15, 2014
% v2.1: Mar. 30, 2014
% v3.2: Aug. 6, 2014
% v4.0: Sept. 13, 2014
% v4.3: Feb. 14, 2015
% v4.4: Jun. 17, 2015
% *************************************************************************


spect = @(x) totalsim(molfreq, obsfreq, aco, agu, eup, lc_idx,...
                      chwidth, beamfilling, x, comp_num);

% Residual function to be minimized.
fitfunc = @(x) norm((sum(spect(x),2) + update - obsint).*nonflagidx);
