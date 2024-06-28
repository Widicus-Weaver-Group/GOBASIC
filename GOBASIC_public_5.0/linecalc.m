% *************************************************************************
% Generate the vectors needed for tau and T_a calculation from catalogs
% Developer: Luyao Zou
% v3.2: Sept. 17, 2014
% v4.3: Feb. 15, 2014

% Input: molcat - molecular catalog, matrix
%        part - molecular partition function file, matrix
% Output: aco - partition function coefficient, double; Q(T) = aco*T^(3/2)
%         agu - product of Einstein coefficient and upper state degeneracy
%               A*g_u, vector
%         eup - upper state energy in J, vector
% *************************************************************************

function [agu, eup] = linecalc(molcat,part)

% Fundamental Constants
k = 1.3806488e-23;  % J/K
c = 29979.2458;     % 10^6 cm
h = 6.62606957e-28; % uJs

% Find the row of T=300
[row, col] = find(part==300);
% Calculate A*gu (Hz) and E_up (J)
int = 10.^molcat(:,3);
% Convert lower state frequency from cm-1 to energy J
elow = h.*c.*molcat(:,5);
eup = elow + h.*molcat(:,1);
agu = 2.7964e-22.*int.*molcat(:,1).^2.*part(row,col-1)...
       ./(exp(-elow/(k*300)) - exp(-eup/(k*300)));
