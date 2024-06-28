% *************************************************************************
% Description: Define varibles and paths for a fitting job
% Author: Mary L. Radhuber
% Date: April 24, 2013
% Developer: Luyao Zou
% v2.0: Dec. 25, 2013
% v2.1: Mar. 28, 2014
% v3.2: Aug. 7, 2014
% v4.1: Sept. 13, 2014
% v4.3: Feb. 14, 2015 HAPPY VALENTINE'S DAY
% v4.4: Jun. 18, 2015
% v4.5: Dec. 1, 2015
% v5.0: Feb. 3, 2017
% *************************************************************************

% >>>>>>>>>>>> OBSERVATIONAL DATA SECTION >>>>>>>>>>>>
%
% Set number of components to fit
comp_num = 2;
%
% Path of observational x,y data
obspath = '../ObsFinal_20141002/W75N.csv';
% Aperture efficiency of the telescope 
% % Except it's actually main beam efficiency, see analysisscript.m line 44
% % main beam Temp = antenna Temp / main beam efficiency
% % - Catherine Walker
apeff = 1;
%
% This 'range' is used to generate 'jlabel'. The value is determined
% by the user, which stands for how many datapoints on each side of the 
% molecular line center are going to be used to calculate gaussian line
% profile. Usually 10*fwhm/(sampling interval) is enough.
% DO NOT CHANGE this value! Changing it will cause error in including 
% previous fitting results.
range = 200;
%%
% Source sizes for each component. If 0, no beam-filling factor is accounted.
% If the source is not circular, please enter the HPFD of the equivilant
% circular Gaussian profile.
% Unit in arcsec. Numbers of source sizes should match 'comp_num'
sourcesize = zeros(comp_num,1);
% Dish diameter of the telescope in meter
dish = 10.4;
%
% Frequency ranges to be flagged. Enter frequency window in MHz row by row.
% Leave blank as "flagfreq = [];" if no flagging is needed.
flagfreq = [% flag 
            %230400,230900;
            %250000,250100;
           ];
%
% Propogate uncertantity from (true) fit,
% (false) noise level before deconvolution
ta_err_from_fit = true; % Default is true
% noise level of antenna temperature (K)
ta_err_deconv = .03; % discarded if ta_err_from_fit = true

% >>>>>>>>>>>> CATALOG SECTION >>>>>>>>>>>>
%
% 'molpath' define catalog paths. Must match 'comp_num'
% 'parpath' define partition function data paths. Must match 'comp_num'
% 'label' define labels for each components for output. Must match 'comp_num'
% 'label' can be left blank if one is too lazy to type in any names.

molpath = {
           '../CatalogFile/CH3OH.txt';
           '../CatalogFile/CH3OH.txt';
           %'more files';
           };

parpath = {
           '../CatalogFile/CH3OHpar.par';
           '../CatalogFile/CH3OHpar.par';
           %'more files';
           };

parco = {
         '../CatalogFile/CH3OHpar.co';
         '../CatalogFile/CH3OHpar.co';
         %'more files';
        };

label = {
         'CH3OH_1';
         'CH3OH_2';
         %'more names';
         };
%
% The transition frequency uncertanty cut-off in MHz from the catalog.
% Frequencies with uncertainty exeeding this cut-off will not be imported.
% 0 if no cut-off desired.
freq_uc_cutoff = 0;
% The log intensity cut-off from the catalog. Transistions weaker than this
% cut-off will not be imported. Should be a negative number.
% 0 if no cut-off desired.
int_cutoff = 0;
% The lower state energy cut-off from the catalog. Lower state energy larger
% than the cut-off will not be imported. Unit in wavenumbers.
% 0 if no cut-off specified.
elow_cutoff = 0;

% >>>>>>>>>>>> FIT PARAMETER SECTION  >>>>>>>>>>>> 
%
% Options for 'patternsearch' algorithm. Customize when necessary.
options = psoptimset('Display','iter','MaxIter',100000,'MaxFunEval',100000*comp_num*5,...
    'ScaleMesh','on','UseParallel','always','CompletePoll','On',...
    'TolMesh',1e-10,'TolX',1e-10,'TolFun',1e-10,'InitialMeshSize',1);

% Boundaries and initial guess of fitted parameters. Must match 'comp_num'
% The four paramters are:
% log(NT) (log10 cm^(-2)), FWHM (km/s), 0.01*Temp (100K), Shift (km/s)
low = repmat([10;2;0;-30],1,comp_num);

up = repmat([20;30;10;30],1,comp_num);

init = repmat([15;15;5;0],1,comp_num);

% The fraction used in calculating uncertainty
tol = 1e-3;

% >>>>>>>>>>>> OUTPUT & PREVIOUS FITS SECTION  >>>>>>>>>>>> 
%
% Define output path and file name. 
outputpath = '../Fitv5.0/HifiBand2/W75N/'; % Don't forget the last '/'
outputname = 'CH3OH_2comp_trial1';
%
% Include previous fits. Only compatible with v4.0 or later versions
pfitpath = {%'../Fitv5.0/HifiBand2/W75N/W75N_FINAL_IX-MOL_Fit.csv';
            %etc;
           };
