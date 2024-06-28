% *************************************************************************
% Description: Generate indices of molecular line frequencies 
%              in the obsfreq vector
% Author: Luyao Zou
% Date: May 30, 2014
% *************************************************************************

function jlabel = jlabel_gen(molfreq,obsfreq,range)
indicies = zeros(1,length(molfreq));
for n = 1:length(molfreq)
    dfreq = abs(obsfreq - molfreq(n));
    index = find(dfreq == min(dfreq));
    indicies(1,n) = index(1);
end 
range_min = -idivide(int16(range),2,'ceil');
range_max = idivide(int16(range),2,'ceil');
range_matrix = repmat((range_min:range_max)',1,length(molfreq));
jlabel = repmat(indicies,range_max-range_min+1,1)-double(range_matrix);
