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
% v5.0: Feb. 3, 2017
% *************************************************************************

% Input: molfreq - molecular transition frequency, vector
%        obsfreq - observational frequency, vector
%        aco - coeff for the partition function
%        agu - einstein coeff * upper state degeneracy
%        eup - upper state energy
%        lc_idx - Matlab index label of line centers for spectrum synthesis
%        chwidth - used to calculate shifted line center indicies
%        x - input parameter vector, contains:
%          x(1) lognt - 10 base logarithm column density in cm^-2
%          x(2) fwhm_v - full width half maximum in km/s
%          x(3) temp*0.01 - in K
%          x(4) shift_v - shift in km/s
% Output: tau - tau(nu) vector over the whole frequency range, sparse
% see 'GOBASICAlgorithm.pdf' for more info

function taumat = tausim(molfreq,obsfreq,aco,agu,eup,lc_idx,chwidth,x)

h = 6.62606957e-34;
k = 1.3806488e-23;
c = 2.99792458e8;
comp_num = size(x,2);

taumat = zeros(length(obsfreq),comp_num);
% sum up tau vector for all molecular components
for n = 1:comp_num
    % prepare variables
    nu = molfreq{n}*1e6;
    nt = 10^(x(1,n)+4);
    fwhm = x(2,n)*molfreq{n}/c*1e3;
    temp = 100*x(3,n);
    shift = x(4,n)*molfreq{n}/c*1e3;
    linecenter = molfreq{n} - shift;
    a = aco{n};
    q = a(1)*temp^(a(2))*(a(3)+exp(a(4)/temp));
    % number of points covered for each Gaussian function
    range_len = size(lc_idx{n}, 1);
    % tau at line center
    tau_c = (sqrt(log(2)/pi^3)/4).*(c.^2./(nu.^2)).*nt./q.*agu{n}.*...
             exp(-eup{n}./(k.*temp)).*(exp(h.*nu./(k.*temp))-1);
    linecenter = repmat(linecenter',range_len,1);
    fwhm = repmat(fwhm',range_len,1);
    shift_idx = repmat(round(shift./chwidth)',range_len,1);
    jlabel = lc_idx{n} - shift_idx;
    ilabel = 1:length(molfreq{n});
    ilabel = repmat(ilabel,range_len,1);
    % construct tau matrix with Gaussian line-profile
    tau = repmat(tau_c',range_len,1)./fwhm.*...
          exp(-4*log(2).*(obsfreq(jlabel)-linecenter).^2./fwhm.^2);
    % reshape tau matrix into a sparse matrix
    % tracing all frequency points in obsfreq
    taulines = sparse(reshape(ilabel,1,[]), reshape(jlabel,1,[]),...
                      reshape(tau,1,[]), length(molfreq{n}), length(obsfreq));
    % sum tau into a vector at the length of obsfreq
    % then store in the n-th column of taumat
    taumat(:,n) = full(sum(taulines,1)');
end
