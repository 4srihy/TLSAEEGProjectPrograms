%fig3
clear;
%get data

load('068TP_F1-200717-GAV_0002_elec30.mat');
elecNum = 1;   elecTrial = 122; %trialNumber
data1 = squeeze(eegData(elecNum,elecTrial,:));
blPos = timeVals>=-0.5 & timeVals<=0;
data = data1(blPos);

t=timeVals(blPos);
freqHFD = [] ;
HFDAll = [];

timeValsTF = [];
freqValsTF = [];
hfdTF = [];

freqRanges = {[1 48] [52 98] [102 148] [152 198] [202 248] [252 298]};
numFreq = length(freqRanges);
colorNames = jet(2*numFreq);
colorNames = flip(colorNames);
colorNames(1:numFreq,:) = colorNames(2:2:2*numFreq,:);
Fs = 2500;
kmax = [10 5];

params.tapers   = [1 1];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 300];
params.trialave = 0;
[SpecPower, freqVals] = mtspectrumc(data',params);

for iFreq = 1:numFreq
    b{iFreq} = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',freqRanges{iFreq}(1),'PassbandFrequency2',freqRanges{iFreq}(2),'PassbandRipple',0.2,'SampleRate',Fs);
    dataFilt{iFreq} = filtfilt(b{iFreq},data')';
    SpecPowerFilt{iFreq} = mtspectrumc(dataFilt{iFreq}',params);
    freqHFD(iFreq) = (freqRanges{iFreq}(1)+freqRanges{iFreq}(2))/2;
    [~,~,~,lnK(iFreq,:),lnL(iFreq,:)] = HigFracDimV2(dataFilt{iFreq},kmax(1),0.05,1);
    hfd(iFreq) = HigFracDimV2(dataFilt{iFreq},kmax(2),0.05,1);
end
lnL0 = lnL - lnL(:,1);
%% time series
figure;

hPlot = getPlotHandles(2,3,[0.07 0.1 0.88 0.86],0.05,0.11);
hPlotFreqTime = getPlotHandles(numFreq/2,2,[0.3800 0.5750 0.2600 0.3650],0.02,0.03);

axes(hPlot(1));
plot(t,data,'m');
xlabel('Time (s)');
ylabel('Voltage (\muV)');
set(hPlot(1),'FontWeight','Bold','FontSize',14);
h=gca;
h.YLim = [h.YLim(1)-30 h.YLim(2)+10];
hPos = h.Position;

hInset = axes('position',[hPos(1)+0.12 hPos(2)+0.03 0.13 0.09]);
plot(timeVals,data1-mean(data1),'Color',[0.6 0.6 0.6]); hold on;
xlim([-0.6 1]);
set(hInset,'FontWeight','Bold','FontSize',10);
rectangle('Position',[-0.5 hInset.YLim(1)+5 0.5 hInset.YLim(2)-hInset.YLim(1)-10],'EdgeColor','m','LineWidth',2);

% timeSeries filtered
for iFreq = 1:numFreq
    plot(hPlotFreqTime(iFreq),t,dataFilt{iFreq}','color',colorNames(iFreq,:),'LineWidth',1); hold on;
    set(hPlotFreqTime(iFreq),'FontSize',12,'FontWeight','bold','TickLength',[0.01 0.01]);
    xlim(hPlotFreqTime(iFreq),[-0.5 0]);
    if iFreq~=1;     ylim(hPlotFreqTime(iFreq),[-5 5]);   end
    title(hPlotFreqTime(iFreq),[num2str(freqRanges{iFreq}(1)) '-' num2str(freqRanges{iFreq}(2))]);

    if iFreq == 2
        ylabel(hPlotFreqTime(iFreq),'Voltage (\mu V)','FontSize',15);
    end

    if ~ismember(iFreq,[numFreq/2 numFreq])
        set(hPlotFreqTime(iFreq),'xTickLabel',[]);
    end

    if iFreq == 3
        axes(hPlotFreqTime(iFreq))
        text(-0.05,-10,'Time (s)','FontSize',15,'FontWeight','bold');
    end

    axes(hPlot(5));
    plot(lnK(iFreq,:),lnL0(iFreq,:),'--','color',colorNames(iFreq,:),'LineWidth',2) ;
    plot(lnK(iFreq,1:kmax(2)),lnL0(iFreq,1:kmax(2)),'color',colorNames(iFreq,:),'LineWidth',3); hold on;% includes dc correction
    legendData{2*iFreq-1} = num2str(round(hfd(iFreq),3));
    legendData{2*iFreq} = {''};
end
legend(hPlot(5),string(legendData),'Location','Southeast');
xlabel(hPlot(5),'ln (1/k)');
ylabel(hPlot(5),'ln (L)');
set(hPlot(5),'FontWeight','Bold','FontSize',14,'TickLength',[0.01 0.01]);

plot(hPlot(2),freqHFD,hfd,'m','LineWidth',2);
xlabel(hPlot(2),'Frequency (Hz)');
ylabel(hPlot(2),'HFD');
set(hPlot(2),'FontWeight','Bold','FontSize',14,'TickLength',[0.01 0.01]);


hPlotTF = getPlotHandles(1,2,[0.3800 0.1000 0.5700 0.3750],0.07,0,0);

%done in the paper
%genTFPlotForHFD(eegData([ 24 26 29 30 31 57 58 61  62 63],:,:),timeVals,1:4,'unipolar',100,8,300,1,'both',hPlotTF); 

%on github..since data shown for one electrode only
genTFPlotForHFD(eegData(1,:,:),timeVals,1:4,'unipolar',100,8,300,1,'both',hPlotTF);

C = findall(gcf,'type','ColorBar');
C(2).Label.String = 'HFD';
C(1).Label.String = '\DeltaHFD';

d1 = colormap("jet");
d2 = colormap(flipud(d1));

for i=1:2
    axes(hPlotTF(i));
    xlabel('Time (s)')
    ylabel('Frequency (Hz)');
    title('');
    set(hPlotTF(i),'FontWeight','Bold','FontSize',14);
end
ylim([0 100]);
a = annotation('textbox',[.0088 0.9415 0.0228 0.0579],'String','(A)','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.3275 0.9415 0.0228 0.0579],'String','(B)','FontSize',16,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.6375 0.9415 0.0228 0.0579],'String','(C)','FontSize',16,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[.0088 0.48 0.0228 0.0579],'String','(D)','FontSize',16,'EdgeColor','none','FontWeight','bold');
e = annotation('textbox',[.3275 0.48 0.0228 0.0579],'String','(E)','FontSize',16,'EdgeColor','none','FontWeight','bold');
f = annotation('textbox',[.6375 0.48 0.0228 0.0579],'String','(F)','FontSize',16,'EdgeColor','none','FontWeight','bold');





