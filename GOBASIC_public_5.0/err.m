% *************************************************************************
% Description: Returns error estimation of the four fitting parameters
% Author: Luyao Zou
% v3.2: Aug. 25, 2014
% v4.3: Feb. 21, 2014
% v4.4: Jun. 13, 2015
% v4.5: Dec. 1, 2015
% *************************************************************************

% Input: x - Parameter vector [lognt, dv, centtemp, shift]
%        tol - percentage tolerance \delta, double
%        ta_err - uncertainty of antenna temperature, double
%        nonflagidx - index for nonflagged values, boolean vector
%        residual - vector
%        spect - spectrum generating function 
%        
% Output: err - associated uncertainty of four parameters
% 
function err = err(x,tol,residual,ta_err,molfreq,obsfreq,aco,agu,...
                   eup,lc_idx,chwidth,beamfilling)

h = 6.62606957e-34;
k = 1.3806488e-23;
coeff = h*1e6/k;

% --- get matrix dimensions
% number of components
comp_num = size(x,2);
% number of variables
numvar = size(x,1) * comp_num;
% length of spectrum
spectlen = length(obsfreq);

% The estimation is time-consuming, since the time complexity is O(n^3)
% [O(n^2) matrix elements, each element calculate O(n) tau values
% however, if we cache the tau for each parameter +/-, time complexity
% can be reduced to O(n^2)
% get spectrum for each component at x
spectcenter = totalsim(molfreq,obsfreq,aco,agu,eup,lc_idx,chwidth,...
                       beamfilling,x,comp_num);

% replicate spectcenter to build up caches
spectcache_p = repmat(spectcenter, 1, 1, size(x,1));
spectcache_m = repmat(spectcenter, 1, 1, size(x,1));
xplus = repmat(x, 1, 1, numvar);
xminus = repmat(x, 1, 1, numvar);
% calculate tau+/tau- at xplus, xminus for each parameter
for n = 1:numvar
    % get the index of n-th parameter in x
    [row, col] = ind2sub(size(x), n);
    xplus(row,col,n) = x(row,col) + tol/2;
    xminus(row,col,n) = x(row,col) - tol/2;
    spectcache_p(:,col,row) = totalsim({molfreq{col}},obsfreq,{aco{col}},...
                              {agu{col}},{eup{col}},{lc_idx{col}},chwidth,...
                              beamfilling(:,col),xplus(:,col,n),1);
    spectcache_m(:,col,row) = totalsim({molfreq{col}},obsfreq,{aco{col}},...
                              {agu{col}},{eup{col}},{lc_idx{col}},chwidth,...
                              beamfilling(:,col),xminus(:,col,n),1);
end

% --- calculate full spectrum at input parameter x
sc = sum(spectcenter, 2);

% calculate first derivative d1
sp = zeros(spectlen, numvar);
sm = zeros(spectlen, numvar);
for n = 1:numvar
    % get the index of the n-th parameter in x
    [row, col] = ind2sub(size(x), n);
    % replace spectcenter with the updated component
    splus = spectcenter;
    sminus = spectcenter;
    splus(:,col) = spectcache_p(:,col,row);
    sminus(:,col) = spectcache_m(:,col,row);
    % cache sp/sm
    sp(:,n) = sum(splus, 2);
    sm(:,n) = sum(sminus, 2);
end
clearvars splus sminus;
d1 = (sp - sm).*repmat(not(not(residual)),1,numvar)./tol;

% calculate second derivative matrix d2
% d2 matrix is symmetric so we only need to calculate half of it
d2 = zeros(numvar);
% get diagonal term of d2
for n = 1:numvar
    d2(n,n) = 4*sum(residual.*(sp(:,n)+sm(:,n)-2.*sc))/(tol^2);
end
clearvars sp sm;

% get off-diagonal term, where 2 parameters are modified simultaneously
% the 4 points:  x(pp) + x(mm) - x(pm) - x(mp)

for m = 1:(numvar-1)
    % get the index of m-th parameter in x
    [rowm, colm] = ind2sub(size(x), m);
    for n = (m+1):numvar
        % get the index of n-th parameter in x
        [rown, coln] = ind2sub(size(x), n);
        % prepare 4 x matrices for 4 x-tol points
        spectpp = spectcenter;
        spectmm = spectcenter;
        spectpm = spectcenter;
        spectmp = spectcenter;
        xpp = x;
        xmm = x;
        xpm = x;
        xmp = x;
        xpp(rowm,colm) = x(rowm,colm) + tol/2;
        xpp(rown,coln) = x(rown,coln) + tol/2;
        xmm(rowm,colm) = x(rowm,colm) - tol/2;
        xmm(rown,coln) = x(rown,coln) - tol/2;
        xpm(rowm,colm) = x(rowm,colm) + tol/2;
        xpm(rown,coln) = x(rown,coln) - tol/2;
        xmp(rowm,colm) = x(rowm,colm) - tol/2;
        xmp(rown,coln) = x(rown,coln) + tol/2;
        % if colm == coln, i.e. from same component
        % x and tau need to be recalculated
        if colm == coln
            spectpp(:,colm) = totalsim({molfreq{colm}},obsfreq,{aco{colm}},...
                              {agu{colm}},{eup{colm}},{lc_idx{colm}},...
                              chwidth,beamfilling(:,colm),xpp(:,colm),1);
            spectmm(:,colm) = totalsim({molfreq{colm}},obsfreq,{aco{colm}},...
                              {agu{colm}},{eup{colm}},{lc_idx{colm}},...
                              chwidth,beamfilling(:,colm),xmm(:,colm),1);
            spectpm(:,colm) = totalsim({molfreq{colm}},obsfreq,{aco{colm}},...
                              {agu{colm}},{eup{colm}},{lc_idx{colm}},...
                              chwidth,beamfilling(:,colm),xpm(:,colm),1);
            spectmp(:,colm) = totalsim({molfreq{colm}},obsfreq,{aco{colm}},...
                              {agu{colm}},{eup{colm}},{lc_idx{colm}},...
                              chwidth,beamfilling(:,colm),xmp(:,colm),1);
        % else if colm ~= coln, i.e. from different components
        % replace spect from spectcache
        else
            spectpp(:,colm) = spectcache_p(:,colm,rowm);
            spectpp(:,coln) = spectcache_p(:,coln,rown);
            spectmm(:,colm) = spectcache_m(:,colm,rowm);
            spectmm(:,coln) = spectcache_m(:,coln,rown);
            spectpm(:,colm) = spectcache_p(:,colm,rowm);
            spectpm(:,coln) = spectcache_m(:,coln,rown);
            spectmp(:,colm) = spectcache_m(:,colm,rowm);
            spectmp(:,coln) = spectcache_p(:,coln,rown);
        end
        % generate spectra
        spp = sum(spectpp,2);
        smm = sum(spectmm,2);
        spm = sum(spectpm,2);
        smp = sum(spectmp,2);
        % calculate off-diagonal matrix elements
        d2(m,n) = sum(residual.*(spp + smm - spm - smp))/(tol^2);
    end
end
% construct full d2 matrix: triangle + t.T - diagonal(double counting)
d2 = d2 + d2' - diag(diag(d2));
clearvars spp smm spm smp taupp taumm taupm taump;

% construct curvature matrix
c = zeros(numvar);
for m = 1:numvar
    for n = 1:numvar
        c(m,n) = (sum(d1(:,m).*d1(:,n)) - d2(m,n))/(ta_err^2);
    end
end
% invert c
epsilon = inv(c);

% uncertainty vector
u = zeros(numvar,1);
for n = 1:numvar
    u(n) = norm(sum(d1.*repmat(epsilon(n,:),spectlen,1),2))/ta_err;
end

% reshape uncertainty to match 'x'
err = reshape(u, size(x,1), comp_num);
