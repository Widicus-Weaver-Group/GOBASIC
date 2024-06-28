% *************************************************************************
% Description: Display results on Matlab terminal and make plots
% Developer: Luyao Zou
% v4.1: Sept. 15, 2014
% v4.4: Jun. 17, 2015
% v4.5: Dec. 1, 2015
% *************************************************************************

% check if label is empty. 
if isempty(label)
    label = cell(comp_num,1);
    for n=1:comp_num
        label{n} = strcat('fit',int2str(n));
    end
end

% Display fitted parameters for each molecular components
for n=1:comp_num
    fprintf('------------- %s --------------\n', label{n})
    fprintf(' NT    = %.4e +/- %.1e cm^-2\n', results(9:10,n))
    fprintf(' FWHM  = %.3f +/- %.3f km/s\n', results(3:4,n))
    fprintf(' Temp  = %.2f +/- %.2f K\n', results(5:6,n))
    fprintf(' shift = %.3f +/- %.3f km/s\n\n', results(7:8,n))
end

% Create label
xlabel({'Frequency (MHz)'});
ylabel({'Intensity (K)'});

% Generate random color specifier
randcolor = rand(comp_num+1,3);

% Make plot
legendstr = cell(1,comp_num+3);
legendstr{1} = 'obs int';
legendstr{2} = 'total fit';
legendstr{3} = 'previous fit';
plot(obsfreq,obsint,'Color',[0.6,0.6,0.6],'linewidth',1.5)
hold on
plot(obsfreq,totalfit,'c-','linewidth',1.5)
plot(obsfreq,update,'Color',randcolor(1,:))
for n=1:comp_num
    legendstr{n+3} = label{n};
    plot(obsfreq,fits(:,n),'Color',randcolor(n+1,:))
end
legend(legendstr,'Location','northeast')

