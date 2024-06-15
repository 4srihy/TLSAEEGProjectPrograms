% fig1 of paper(v5 at original file)
% adding HFD computation for artificial data like straight line, sin, noisy
% sin, gaussian, pink noise, wfbm
% to add p values in violin plots
figure;

% Load the data for human single trial for elec 30
load('068TP_F1-200717-GAV_0002_elec30.mat');
elecNum = 1;   elecTrial = 122; %trialNumber
dataAllTrials = squeeze(eegData(elecNum,:,:));
data1 = squeeze(eegData(elecNum,elecTrial,:));

colorNamesBl = [1 0 1; 0 1 1; 0.5 0.8 0.7];
lineColor = [0.5 0.5 0.5];%[1 0.9 0.5];

numPlots = 6;
colorNames = turbo(2*numPlots);
colorNames = flip(colorNames);
colorNames(1:numPlots,:) = colorNames(2:2:2*numPlots,:);

blRange = [-0.5 0];
stRange = [0.25 0.75];
blPos = timeVals>=blRange(1) & timeVals<=blRange(2);
stPos = timeVals>=stRange(1) & timeVals<=stRange(2);

freqRange = [1 90];
Fs = 2500;
kmax = [10 30];


b = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',freqRange(1),'PassbandFrequency2',freqRange(2),'PassbandRipple',0.2,'SampleRate',Fs);
dataFiltInd = filtfilt(b,data1')';
dataFiltAll = filtfilt(b,dataAllTrials')';
dataBl = dataFiltInd(blPos);
dataSt = dataFiltInd(stPos);

[hfdBl,~,~,lnKBl,lnLBl] = HigFracDimV2(dataBl,kmax(1),0.05,1);
[hfdSt,~,~,lnKSt,lnLSt] = HigFracDimV2(dataSt,kmax(1),0.05,1);

[~,~,~,lnKBl2,lnLBl2] = HigFracDimV2(dataBl,kmax(2),0.05,1);
[~,~,~,lnKSt2,lnLSt2] = HigFracDimV2(dataSt,kmax(2),0.05,1);

dcCorrecLnL = mean([lnLBl(1) lnLSt(1)]);
dcCorrecLnL2 = mean([lnLBl2(1) lnLSt2(1)]);

for iTrial = 1:size(dataAllTrials,1)
    hfdBlAll(iTrial) = HigFracDimV2(dataFiltAll(iTrial,blPos),kmax(1),0.05,1);
    hfdStAll(iTrial) = HigFracDimV2(dataFiltAll(iTrial,stPos),kmax(1),0.05,1);
    hurstBl(iTrial) = genhurst(dataFiltAll(iTrial,blPos));
    hurstSt(iTrial) = genhurst(dataFiltAll(iTrial,stPos));
end

erp = mean(dataAllTrials-mean(dataAllTrials,2),1);

hPlot = getPlotHandles(1,4,[0.07 0.59 0.9 0.3750],0.05,0.11,0);

%for Synthetic signals
hPlot2 = getPlotHandles(3,2,[0.0700 0.59 0.1875 0.3750],0.02, 0.03,0);
totIter = 10; % no. of iterations for synthetic noise
nPts = 20000;
ptsUsed = 1000; % no. of data points to use in analysis
peakfreq = 10;
amplitude = 0.5;
H = 0.3; % only if series is FBM; H>.5 deterministic
m = 0.5;
slope = 1.5;
t = (1:nPts)/2500;

% for slope calculation for fBm : works only if fooof is installed
% settings = struct();
% settings.max_n_peaks = 0;
% fres = round(1/(max(t(ptsUsed)+1/Fs)));
% f = 0:fres:(ptsUsed-1)*fres;
% freqRangeSlope = [1 f(ptsUsed/2)];

for iter = 1:totIter
    x{1}{iter} = m.*t; %straight line

    x{2}{iter} = amplitude*sin(2*pi*peakfreq*t); % pure sin

    x{3}{iter} = amplitude*sin(2*pi*peakfreq*t) + 0.4*amplitude*rand(1,nPts); %noisy sinusoid

    freqRes = 2500/nPts;
    freq = freqRes:freqRes:1250;
    power = 10^4./freq.^(slope);
    fftDataAmp([1:nPts/2 nPts:-1:nPts/2+1]) = [sqrt(power) sqrt(power)];
    phaseData0 = -3.14 + 6.28*rand(1,nPts/2);
    phaseData([1:nPts/2 nPts:-1:nPts/2+1]) = [phaseData0 -phaseData0];
    fftDataFull = fftDataAmp.*(exp(1i.*phaseData));

    %do ifft
    x{4}{iter} = ifft(fftDataFull,'symmetric');% colored noise

    x{5}{iter} = wfbm(H,nPts);

    x{6}{iter} = rand(1,nPts);
end

syntheticList = [{'Straight Line'},{'Pure Sinusoid'},{'Noisy Sinusoid'},{'Colored Noise'},{'FBM (H=0.3)'},{'Gaussian Noise'}];

for i=1:numPlots
    for iter=1:totIter
        [hfd0(i,iter),~,~,~,logL0{i}(iter,:)] = HigFracDimV2(x{i}{iter}(1:ptsUsed),10,0.05,0,0);
        dcCorr0(iter) = logL0{i}(iter,1);

        if i==4
            hurst(iter) = genhurst(x{i}{iter}(1:ptsUsed));
            %     elseif i==5
            %         SpecPower = (abs(fft(x{i}{iter}(1:ptsUsed)))).^2;
            %         fooof_ind = fooof(f,SpecPower,freqRangeSlope,settings,true);
            %         exponent(iter) = fooof_get_params(fooof_ind,mean(SpecPower,1)',settings);
        end
    end

    logL1{i} = mean(logL0{i}-(logL0{i}(:,1)-mean(dcCorr0)));
    stdL1{i} = std(logL0{i}-(logL0{i}(:,1)-mean(dcCorr0)));
    dcCorr(i) = logL1{i}(1);
end
hfd = mean(hfd0,2);
dcCorrMean = mean(dcCorr);

for i=1:numPlots
    plot(hPlot2(i),t(1:1000),x{i}{1}(1:1000),'color',colorNames(i,:),'LineWidth',1); hold on;
    set(hPlot2(i),'FontSize',12,'FontWeight','bold','TickLength',[0.01 0.01]);
    title(hPlot2(i),syntheticList{i});
    if ~ismember(i,[numPlots/2 numPlots])
        set(hPlot2(i),'xTickLabel',[]);
    end
    if i == 2
        ylabel(hPlot2(i),'Voltage (\muV)','FontSize',15);
    end

    if i== 3
        axes(hPlot2(i))
        text(0.37,-1.2,'Time (s)','FontSize',15,'FontWeight','bold');
    end

    axes(hPlot(2));
    xs = -log(1:10);
    y = logL1{i}-(logL1{i}(1)-dcCorrMean);
    plot(-log(1:10),logL1{i}-(logL1{i}(1)-dcCorrMean),'color',colorNames(i,:),'LineWidth',1.5); hold on;
    patch([xs';flipud(xs')],[y'-stdL1{i}';flipud(y'+stdL1{i}')],colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
    legendData{2*i-1} = num2str(round(hfd(i),3));
    legendData{2*i} = {''};
end
legend(hPlot(2),string(legendData),'Location','Southeast');
xlabel(hPlot(2),'log (1/k)');
ylabel(hPlot(2),'log (L)');
set(hPlot(2),'FontWeight','Bold','FontSize',14,'TickLength',[0.01 0.01]);

axes(hPlot(3));
plot(timeVals,dataFiltInd,'LineWidth',2,'Color',lineColor); hold on;
xlim([-0.6 1]);
h=gca;
makeBox(hPlot(3),blRange,h.YLim,colorNamesBl(1,:),2.5,'--','V');
makeBox(hPlot(3),stRange,h.YLim,colorNamesBl(2,:),2.5,'--','V');
set(hPlot(3),'FontWeight','bold','FontSize',14,'Box','off');
xlabel('Time (s)','FontSize',14,'FontWeight','Bold');
ylabel('Voltage (\muV)','FontSize',14,'FontWeight','Bold');


axes(hPlot(4))
plot(lnKBl,lnLBl-(lnLBl(1)-dcCorrecLnL),'color',colorNamesBl(1,:),'LineWidth',2); hold on;
plot(lnKSt,lnLSt-(lnLSt(1)-dcCorrecLnL),'color',colorNamesBl(2,:),'LineWidth',2);
set(hPlot(4),'FontWeight','bold','FontSize',14,'Box','off');
xlabel('log (1/k)','FontSize',14,'FontWeight','Bold');
ylabel('log (L)','FontSize',14,'FontWeight','Bold');
legend(['BL(' num2str(round(hfdBl,3)) ')'],['ST(' num2str(round(hfdSt,3)) ')'],'Location','northwest');


%% to plot hfd for all subjects
folderName = fullfile(folderSourceString{1},projectName,'higFracDim');
stFlag = 1;
medianFlag=1;
diseaseFlag = 0;
dimToMean = 2;
groupNum = 6;
clear dim0 dimS0 dimDiff parameter hurst0 hurstS0
for i=1:length(subjectNameList)
    for j=1:length(subjectNameList{i})
        fileName = fullfile(folderName,[subjectNameList{i}{j} '_' refType '.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
            for k=1
                dim0{k}{i}(j,1:64) = nan;
                dimS0{k}{i}(j,1:64) = nan;
                hurst0{k}{i}(j,1:64) = nan;
                hurstS0{k}{i}(j,1:64) = nan;
            end

        else
            xj = load(fileName);

            for k=1:length(xj.fracDim.higBL2)
                dim0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higBL2{k}(:,:,1),dimToMean));
                dimS0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higSt2{k}(:,:,1),dimToMean));
                if k==1
                    hurst0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.hurst2,dimToMean));
                    hurstS0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.hurstSt2,dimToMean));
                    sum0{k}{i}(j,:) = dim0{k}{i}(j,:) + hurst0{k}{i}(j,:);
                end

            end
        end
    end
end

parameter{1}{1} = dim0{8};     parameter{1}{2} = dimS0{8};       parameterName{1} = 'HFD';
parameter{2}{1} = hurst0{1};   parameter{2}{2} = hurstS0{1};      parameterName{2} = 'HI';

%% plotting
if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
    [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    groupNameList{length(groupNameList)+1} = 'highPriority';
end

numGroup = length(subjectNameList);
capType = 'actiCap64';

totCondns = length(parameter);
iElecGroup = groupNum;

for iCond = 1:totCondns
    for iSub = 1:2
        if medianFlag
            dataSwarm{iCond}{iSub} = cat(1,nanmedian(parameter{iCond}{iSub}{1}(:,electrodeGroupList{iElecGroup}),2),nanmedian(parameter{iCond}{iSub}{2}(:,electrodeGroupList{iElecGroup}),2));
        else
            dataSwarm{iCond}{iSub} = cat(1,nanmean(parameter{iCond}{iSub}{1}(:,electrodeGroupList{iElecGroup}),2),nanmean(parameter{iCond}{iSub}{2}(:,electrodeGroupList{iElecGroup}),2));
        end
    end
    dataY{iCond+1} = [dataSwarm{iCond}{1}' dataSwarm{iCond}{2}'];
end

dataY0{1}{1} = hfdBlAll;    dataY0{1}{2} = hfdStAll;
for iu = 1:2
    dataY0{2}{iu} = dataSwarm{1}{iu};
end


%%  add violin plots  % legend off is not working
dataY{1} = [hfdBlAll hfdStAll];
dataY{4} = [hurstBl hurstSt];

for iCond = 1:2
    numDataPoints(iCond) = size(dataY{iCond},2)/2;
    dataDiff{iCond} = dataY{iCond}(numDataPoints(iCond)+1:end)-(dataY{iCond}(1:numDataPoints(iCond)));
end


posx = [ 0.25 , 0.73 ];
yLabels = [{'HFD'}, {'HFD'}, {'H'}];

tempx = [{'All Trials'}, {'All Subjects'}];

clear g
for iCond = 1:2

    %data diff
    xDiff = repelem(tempx,numDataPoints(iCond));
    if iCond==1
        violdata{1}.x = xDiff(1:numDataPoints(iCond));
    else
        violdata{1}.x = xDiff(numDataPoints(iCond)+1:2*numDataPoints(iCond));
    end
    violdata{1}.y = dataDiff{iCond};
    violdata{1}.group = violdata{1}.x;

    if iCond~=1
        violdata{2}.x = repelem({'BL','ST'},numDataPoints(iCond));
        violdata{2}.y = dataY{iCond};
        violdata{2}.group = violdata{2}.x;
    end

    for s = 1

        g(s,iCond)=gramm('x',violdata{s}.x,'y',violdata{s}.y,'color',violdata{s}.group);
        g(s,iCond).stat_violin('fill','transparent','normalization','width');
        g(s,iCond).stat_boxplot('width',0.13);

        if s==1
            g(s,iCond).set_color_options('map',colorNamesBl(3,:));
        else
            g(s,iCond).set_color_options('map',colorNamesBl);
        end
        g(s,iCond).axe_property('FontWeight','bold','FontSize',14);

        g(s,iCond).set_layout_options('position',[posx(s,iCond) 0.012 0.26 0.45],'legend_position',[0 0 0 0]);

        g(s,iCond).set_names('x','','y',yLabels{iCond},'Color','');
        g(s,iCond).set_order_options('x',0,'color',0);

    end
end

g.draw();

%% adding violin plots including pairing lines
displaySettings.parametricTest = 0;
displaySettings.showYTicks=1;
displaySettings.showXTicks=1;
displaySettings.commonYLim = 1;
displaySettings.xPositionText =0.8;
displaySettings.textFontSize = 10;
displaySettings.yPositionLine=0.15;
displaySettings.xTickLabels = [{'BL'} {'ST'}];

for iCond = 1:2
    h1 = axes('position',[0.07+0.4756*(iCond-1) 0.066 0.185 0.38]);
    displaySettings.plotAxes = h1;

    switch iCond
        case 1
            displaySettings.setYLim=[1.01 1.2];
        case 2
            displaySettings.setYLim=[1.04 1.18];
    end
    ylabel(h1,'HFD');
    ax=displayViolinPlot_Sri(dataY0{iCond},[{colorNamesBl(1,:)} {colorNamesBl(2,:)}],1,1,0,1,displaySettings);
    set(ax,'FontWeight','bold','FontSize',14);
end

allAxesInFigure = findall(gcf,'type','axes');
%p -value
diffParam = -[0.005 0.005 0.005];
deltaXPos = [0.03 0.11 -0.02];
deltaYPos = [-0.02 -0.0054 -0.001];
for iCond = 1:2
    % p value Here only for equal sets
    numData = size(dataY{iCond},2)/2;
    if medianFlag
        pD(iCond)=signrank(dataY{iCond}(1:numData),dataY{iCond}(numData+1:end));
    else
        %t-test
        [~,pD(iCond)]=ttest(dataY{iCond}(1:numData),dataY{iCond}(numData+1:end)); %from 'ttest
    end

    h =allAxesInFigure(2*(4-iCond));% allAxesInFigure(2*(6-iCond));
    plot(h,[h.XLim(1) h.XLim(2)-0.2],[0 0],'--k','LineWidth',1);
    pPosition = h.YLim(2)-diffParam(iCond)/5;

    if pD(iCond)<0.001
        t =   text(h,0.8,pPosition,num2str(pD(iCond),'%.1e'),'FontSize',12,'FontWeight','bold');
    else
        text(h,0.9,pPosition,num2str(pD(iCond),'%.3f'),'FontSize',12,'FontWeight','bold');
    end
    text(h,deltaXPos(iCond),deltaYPos(iCond),'\Delta','FontSize',14,'FontWeight','bold','Rotation',90);

end

a = annotation('textbox',[.0144 0.9525 0.0228 0.0579],'String','(A)','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.4957 0.9525 0.0228 0.0579],'String','(B)','FontSize',16,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.0144 0.4553 0.0228 0.0579],'String','(C)','FontSize',16,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[.4957 0.4553 0.0228 0.0579],'String','(D)','FontSize',16,'EdgeColor','none','FontWeight','bold');

