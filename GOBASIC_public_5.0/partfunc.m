% *************************************************************************
% Description: Finds an interpolated roational partition function from the .par file
% Author: Mary Radhuber
% Date: June 15, 2012
% Developer: Luyao Zou
% v2.0: Dec. 11, 2013
% v4.0: Sept. 15, 2014
% v4.5: Dec. 1, 2015
% *************************************************************************

function a = partfunc(part)
partition = @(a,x) a(1)*x.^(a(2));
% 3 sets of initial guess are provided, try each out 
try
    guess = [part(1,1)/(part(1,2)^1.5); 1.5];
    a = nlinfit(part(:,2), part(:,1), partition, guess);
catch ME
    % if ^1.5 does not work, try ^2
    try
        guess = [part(1,1)/(part(1,2)^2); 2];
        a = nlinfit(part(:,2), part(:,1), partition, guess);
    catch ME
        % if it still does not work, try ^1
        try
            guess = [part(1,1)/part(1,2); 1];
            a = nlinfit(part(:,2), part(:,1), partition, guess);
        catch ME
            rethrow(ME)
        end
    end
end
