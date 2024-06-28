% *************************************************************************
% Description: Truncate catalog files into fitted frequency range and given
%             cut-off settings.
% Author: Mary L. Radhuber
% Date: April 24, 2013
% Developer: Luyao Zou
% v2.0: Dec. 25, 2013
% v4.0: Sept. 14, 2014
% *************************************************************************

% Input : mol - catalog file
%         freq_min - frequency window minimum
%         freq_max - frequency window maximum
%         fcut - frequency uncertainty cut-off
%         intcut - intensity cut-off
%         elowcut - upper state energy cut-off
% 
% Output: molcal - truncated catalog data

function molcal = trunccat(mol,freq_min,freq_max,fcut,intcut,elowcut)
% cut-off frequency window
idx = mol(:,1)>freq_min & mol(:,1)<freq_max;

% cut-off uncertainties larger than fcut
if(fcut>0)
    idx = mol(:,2)<=fcut & idx;
end

% cut-off intensity smaller than intcut
if(intcut<0)
    idx = mol(:,3)>=intcut & idx;
end

% cut-off lower state energy larger than elowcut
if(elowcut>0)
    idx = mol(:,5)<=elowcut & idx;
end

% pass truncated catalog table
molcal = mol(idx,:);
