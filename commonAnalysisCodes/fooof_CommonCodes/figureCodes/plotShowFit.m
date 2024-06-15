%figure

%load ()
%load the subject data for which you want to plot it 
elecNum = 31; %the electrode no. for which data is plotted
figure;
h1 = subplot('Position',[0.07 0.1 0.43 0.88]);

for jFreq = 1:length(freq_range)
    fooof_plot(fooof_results{elecNum}(jFreq),true);
end
xlim([4 800]);
legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit')
set(h1,'LineWidth',1,'FontSize',10,'FontWeight','bold');

h2 = subplot('Position',[0.54 0.1 0.43 0.88]);
plot(centreFreq,exponent(32,:),'b',"LineWidth",3)
xlim([30 800]);
set(h2,'XScale','log','LineWidth',1,'FontSize',10,'FontWeight','bold','box','off');

ylabel(h1,'log_{10} (Power (\muV)^2)','FontWeight','bold','FontSize',14);
ylabel(h2,'Exponent','FontWeight','bold','FontSize',14);
xlabel(h1,'Frequency (Hz)','FontWeight','bold','FontSize',14);
xlabel(h2,'Frequency (Hz)','FontWeight','bold','FontSize',14);