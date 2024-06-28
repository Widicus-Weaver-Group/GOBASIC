% *************************************************************************
% Description: Import files and calibrate data
% Author: Mary Radhuber
% Date: May 1, 2013
% Developer: Luyao Zou
% v2.0: Dec. 25, 2013
% v3.0: May 30, 2014
% v3.2: Aug. 6, 2014
% v4.1: Sept. 13, 2014
% v4.3: Feb. 15, 2015
% v4.4: Jun. 13, 2015
% v5.0: Feb. 3, 2017
% *************************************************************************

% Check possible errors
if size(sourcesize,2) ~= comp_num
    error('GOBASIC:CompNumMismatch','Number of source sizes does not match number of components')
end
if length(molpath) ~= comp_num
    error('GOBASIC:CompNumMismatch','Number of catalog files does not match number of components')
end
if length(parpath) ~= comp_num
    error('GOBASIC:CompNumMismatch','Number of partition fucntion files does not match number of components')
end
if ~isempty(label) && length(label)~= comp_num
    error('GOBASIC:CompNumMismatch','Number of labels does not match number of components')
end
if (size(low,2)~=comp_num || size(up,2)~=comp_num || size(init,2)~= comp_num)
    error('GOBASIC:CompNumMismatch','Dimension of parameter matrix does not match number of components')
end

% Import observational data, flagging and calibrate it by apeture efficiency
data = importdata(obspath);
% Retrieve the minimum and maximum frequency of the data
freq_min = data(1,1);
freq_max = data(end,1);
% Extend the data frequency a little bit to prevent possible error.
% This extension allows simulation of transitions near the edge of the
% real data, in which case the gaussian function may goes out of the data
% frequency range. Intensity is filled with 0 so real data stays unmodified.
data = cat(1,[(freq_min-range:freq_min-1)' zeros(range,1)],...
           data,[(freq_max+1:freq_max+range)' zeros(range,1)]);
obsfreq = data(:,1);
if (apeff>0 && apeff<=1)
    obsint = data(:,2)/apeff;
elseif apeff>1
    error('GOBASIC:ValueInvalid','Apperture efficiency cannot exeed 1')
elseif apeff<0
    error('GOBASIC:ValueInvalid','Apperture efficiency cannot be negative')
end

% Used for calculating shifted line-center indicees
% It is the mean value of the adjacent x-data difference, i.e. channel width
chwidth = mean(diff(obsfreq));

% Flag data. The data flagged are omitted during optimization
flagidx = flag(obsfreq,flagfreq);
nonflagidx = not(flagidx);

% The cell array 'fitcomps' stores arrays needed in fitting for all
% molecular components. The size of the cell array changes dynamically
% based on the number of components specified. So theoretically one
% can fit infinit number of components. However, computational cost
% increases exponentially.
%
% These are the cell arrays for the n fitting components
% Cell size changes automatically base on user definition
% aco is the coeff before partition function, Q(T)=aco*T^(3/2)
% agu is the einstein coeff * upper state degeneracy
% molcat is the truncated catalog file
% molfreq is the transition frequencies
% eup is the upper state energy
% jlabel is an Matlab index label used for spectrum synthesis
aco = cell(comp_num,1);
agu = cell(comp_num,1);
molcat = cell(comp_num,1);
molfreq = cell(comp_num,1);
eup = cell(comp_num,1);
lc_idx = cell(comp_num,1);
% Generate arrays for each molecular components
for n=1:comp_num
    tempcat = importdata(molpath{n});
    molcat{n} = trunccat(tempcat,freq_min,freq_max,...
                         freq_uc_cutoff,int_cutoff,elow_cutoff);
    molfreq{n} = molcat{n}(:,1);
    part = importdata(parpath{n});
    %[aco{n}, qrms{n}] = partfunc(part);
    aco{n} = importdata(parco{n});
    [agu{n}, eup{n}] = linecalc(molcat{n},part);
    lc_idx{n} = linecenter(molfreq{n},obsfreq,range);
end
clearvars tempcat part;

% calculate beam filling factor if source size is specified
c = 2.99792458e8;
% get corresponding wavelength at each frequency point
lambda = c./(obsfreq.*1e6);
% beam size at each frequency (arc seconds)
beamsize = (sqrt(-8*log(2)*log(0.1)).*lambda./dish*180*3600/pi^2);
sourcesize = repmat(sourcesize,length(obsfreq),1);
beamsize = repmat(beamsize,1,comp_num);
beamfilling = sourcesize.^2./(sourcesize.^2 + beamsize.^2);
% replace 0 with 1
beamfilling = beamfilling + not(beamfilling);
clearvars sourcesize beamsize dish lambda;

% Load previous fits
% Include previous fits
update = zeros(length(obsint),1);
for n = 1:length(pfitpath)
    pfit = importdata(pfitpath{n});
    update = update + pfit.data(:,end);
end
clearvars pfit;
