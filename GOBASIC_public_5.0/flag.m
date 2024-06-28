% *************************************************************************
% Description: Flag RFI spikes
% Author: Luyao Zou
% v3.2: Aug 6, 2014
% *************************************************************************

% Returns matrix indices of flagged frequencies
function flagidx = flag(obsfreq,flagfreq)
flagidx = false(length(obsfreq),1);
for n=1:size(flagfreq,1)
    flagidx = flagidx + (obsfreq>=flagfreq(n,1) & obsfreq<=flagfreq(n,2));
end
