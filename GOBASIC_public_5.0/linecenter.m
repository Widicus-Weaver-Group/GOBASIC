% *************************************************************************
% Description: Generate indices for frequency centers of molecular lines
%              in the obsfreq vector
% Author: Luyao Zou
% v3.0: May 30, 2014
% *************************************************************************

function linecenter_idx = linecenter(molfreq,obsfreq,range)
indicies = zeros(1,length(molfreq));
for n = 1:length(molfreq)
    dfreq = abs(obsfreq - molfreq(n));
    index = find(dfreq == min(dfreq));
    indicies(1,n) = index(1);
end 
range_min = -idivide(int16(range),2,'ceil');
range_max = idivide(int16(range),2,'ceil');
range_matrix = repmat((range_min:range_max)',1,length(molfreq));
linecenter_idx = repmat(indicies,range_max-range_min+1,1) - double(range_matrix);
