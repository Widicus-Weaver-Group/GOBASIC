% *************************************************************************
% Description: Simulates spectrum
% Author: Mary L. Radhuber
% Date: June 15, 2012
% Developer: Luyao Zou
% v2.0: Dec. 11, 2013
% v3.0: May 30, 2014
% v3.2: Aug. 23, 2014
% v4.0: Sept. 15, 2014
% v4.3: Feb. 27, 2015
% v4.4: Jun. 17, 2015
% v4.5: Dec. 1, 2015
% *************************************************************************

% Input: molfreq - molecular transition frequency, vector
%        obsfreq - observational frequency, vector
%        aco - coeff for the partition function, Q(T)=aco(1)*T^(aco(2))*exp(aco(3)/T)
%        agu - einstein coeff * upper state degeneracy
%        eup - upper state energy
%        lc_idx - Matlab index label of line centers for spectrum synthesis
%        chwidth - used to calculate shifted line center indicies
%        sourcesize - source size in arcsec
%        dish - dish diameter of the telescope in m
%        lognt - 10 base logarithm column density in cm^-2
%        fwhm_v - full width half maximum in km/s
%        lntemp - natural logarithm temperature in K
%        shift_v - shift in km/s
% Output: totalsim - spectrum of a single molecular component, vector
% see 'GOBASICAlgorithm.pdf' for more info

function spectra = totalsim(molfreq, obsfreq, aco, agu, eup, lc_idx,...
                            chwidth, beamfilling, x, comp_num)

h = 6.62606957e-34;
k = 1.3806488e-23;
coeff = h*1e6/k;

tau = tausim(molfreq,obsfreq,aco,agu,eup,lc_idx,chwidth,x);
spectra = coeff.*repmat(obsfreq,1,comp_num).*(1-exp(-tau))./(exp(coeff.*...
          repmat(obsfreq,1,comp_num)./(repmat(100.*x(3,:),length(obsfreq),...
          1)))-1).*beamfilling;
% replace possible NaN values with 0
spectra(isnan(spectra))=0;
