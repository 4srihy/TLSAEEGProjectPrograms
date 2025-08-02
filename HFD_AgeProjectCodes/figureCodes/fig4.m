%figure 4  %31-1-24

%adding HFD variation with Poower

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotFigures; 
diseaseFlag=0;  
bonferroniFlag = 0;   

fooofFlag = 1;
knee_flag = 0;   % for evaluating exponent for subject average psd
powerType = 'BL';
stFlag = 1;
folderName = fullfile(folderSourceString{1},projectName,'higFracDim');

if fooofFlag
    if knee_flag
        str_condition = 'withKnee';
    else
        str_condition = 'withoutKnee';
    end
else
    str_condition = 'withoutFooof';
end

if ~exist('capType','var');  capType = 'Acticap64';  end
showSEMFlag=1;
medianFlag = 1;
OptimR_SQ = 0;
OptimExponent = 0.01;
smoothSigma = 5;
parameter = 'exponent';
cLimParam = [0 4]; % limit of plots
nanValue =0; %the value of parameter
groupNum = 6; %electrode group num to be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%variable parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tapers = 1;%1/2/3[1 2 3];
freqRangeWidths = [100 76 200];
freqRangeList = {[64 140],[230 430]};
freqList = [{'LFR (64-140 Hz)'} {'HFR (230-430 Hz)'}];
freqWidthHfd = 100;
kMax = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the High priority electrodes
%need to be independant

if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
    [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    groupNameList{length(groupNameList)+1} = 'highPriority';
end


%%%%%%%%%%%%%%%% for different groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numGroups = length(subjectNameList);

y=[];
err= [];
if length (folderSourceString)<=3
    totDataType = length(folderSourceString);
else
    totDataType = 1;
    folderSourceStringNew{1} = folderSourceString;
    clear folderSourceString;
    folderSourceString = folderSourceStringNew;
end

dataTypeList{1} = 'Eyes Open'; dataTypeList{2} = 'Eyes Closed';
jDat = 1;

for i = 1:numGroups
    %for continuous freq range
    [freqVals{jDat},centreFreq,PSD{jDat}{i},paramData{jDat}{i}, powerElectrodeGroup{jDat}{i},paramElectrodeGroupData{jDat}{i},paramSubjectElectrodeGroupData{jDat}{i},logPSD{jDat}{i}, logPowerElectrodeGroup{jDat}{i},subjectNum{jDat}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(1),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);

    for j=1:length(subjectNameList{i})
        fileName = fullfile(folderName,[subjectNameList{i}{j} '_' refType '_freqWidth' num2str(freqWidthHfd) '_kmax' num2str(kMax) '.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
        else
            xj = load(fileName);
            for k=1:length(xj.fracDim.higBL2)
                if freqWidthHfd==24
                    zc  = bipolarMean(xj.fracDim.higBL2{k}',0);
                    dim0{i}(j,k,:) = zc(1,:);
                    if stFlag
                        zcs = bipolarMean(xj.fracDim.higSt2{k}',0);
                        dimS0{i}(j,k,:) = zcs(1,:);
                    end
                else
                    dim0{i}(j,k,:) = squeeze(nanmedian(xj.fracDim.higBL2{k}(:,1,:),3));
                    if stFlag;  dimS0{i}(j,k,:) = squeeze(nanmedian(xj.fracDim.higSt2{k}(:,1,:),3));   end
                end
            end
        end
    end
    dim0{i}(dim0{i}==0) = nan;

end

if stFlag
    for i=1:2
        dimS0{i}(dimS0{i}==0) = nan;
        dimDiff{i} = dimS0{i}-dim0{i};
    end
end

if freqWidthHfd == 24
    freqValsHFD = xj.fracDim.centreFreq(1:29);
else
    freqValsHFD = squeeze(xj.fracDim.centreFreq(:,:,1));
end

paramHFD = dim0;

displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
xScaleVar = 'linear'; %Log/Linear
figure;
colormap jet;
colorNames = [1 0 0; 0 0 1; 1 0 1;0.49 0.18 0.56; 0.5 0.5 0.5; 0 1 1 ; 0 1 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%hot(8); %colorNames([1:3,end-2:end],:) = [];
displaySettings.colorNames = colorNames;
color1 = [1 0.75 0.13];    color2=[1 0.3 0]; %add colors for line

displaySettingsDiff = displaySettings;
displaySettingsDiff.colorNames = [0.49 0.18 0.56];

%% figures psd and exponent
axTopo = [];    axParam=[]; axPSD = []; yLimPSD = [];  yLimParam = [];


centreFreqNew = centreFreq;
clear dataPSD dataParam dataTopo dataHFD
for i = 1:numGroups
    if groupNum ==0
        if medianFlag
            dataPSD{i} = squeeze(nanmedian(logPSD{jDat}{i},2));
            dataParam{i} = squeeze(nanmedian(paramData{jDat}{i},2));
            dataHFD{i} = squeeze(nanmedian(paramHFD{i},3));
        else
            dataPSD{i} = squeeze(nanmean(logPSD{jDat}{i},2));
            dataParam{i} = squeeze(nanmean(paramData{jDat}{i},2));
            dataHFD{i} = squeeze(nanmean(paramHFD{i},3));
        end
    else
        dataPSD{i} = logPowerElectrodeGroup{jDat}{i}{groupNum};
        dataParam{i} = paramSubjectElectrodeGroupData{jDat}{i}{groupNum};
        dataHFD{i} = squeeze(nanmedian(paramHFD{i}(:,:,electrodeGroupList{groupNum}),3));
    end
end
clear ylimparam

param{1} = dataHFD;     xVals{1} = freqValsHFD;         paramLabel{1} = 'HFD';      ylimparam{4} = [-0.008 0.008];
param{2} = dataPSD;     xVals{2} = freqVals{jDat};      paramLabel{2} = 'log_{10} (Power (\muV)^2)';
param{3} = dataParam;   xVals{3} = centreFreqNew;       paramLabel{3} = 'Slope';    ylimparam{3} = [0 3];    ylimparam{6} = [-0.4 0.4];

displaySignificanceFlag=1;
bFlag  = bonferroniFlag;

numRows = 5;

for iParam = 1:length(param)

    hPlot(iParam) = subplot(numRows,2,2*(iParam-1)+1);%main
    hPlot(iParam+3) = subplot(numRows,2,2*iParam);%diff

    if isempty(ylimparam{iParam})
        ylim1 = '';
    else
        ylim1 = ylimparam{iParam};
    end

    if isempty(ylimparam{iParam+3})
        ylim2 = '';
    else
        ylim2 = ylimparam{iParam+3};
    end


    displayAndcompareData(hPlot(iParam),param{iParam},xVals{iParam},displaySettings,ylim1,displaySignificanceFlag,medianFlag,'','','',bFlag); xlim([6 300]);

    if iParam==2
        for iNum=1:2;    param{iParam}{iNum} = 10*param{iParam}{iNum};   end  % to display change in power in db
    end
    displayAndcompareData(hPlot(iParam+3),param{iParam},xVals{iParam},displaySettingsDiff,ylim2,displaySignificanceFlag,medianFlag,'','','',bFlag,1); xlim([6 300]); hold on;
    box off;
    plot(hPlot(iParam+3),xVals{iParam},0.*xVals{iParam},'k','LineWidth',1);
    yLim{iParam} = hPlot(iParam).YLim;
    set(hPlot(iParam),'Xscale',xScaleVar,'LineWidth',1,'FontSize',10,'FontWeight','bold');
    ylabel(hPlot(iParam),paramLabel{iParam},'FontWeight','bold','FontSize',12);

    yLim{iParam+3} = hPlot(iParam+3).YLim;
    set(hPlot(iParam+3),'Xscale',xScaleVar,'LineWidth',1,'FontSize',10,'FontWeight','bold');

    if iParam==2
        ylabel(hPlot(iParam+3),'\Delta Power (dB)','FontWeight','bold','FontSize',12);
        set(hPlot(iParam),'XTickLabel',[]);
        set(hPlot(iParam+3),'XTickLabel',[]);
    else
        ylabel(hPlot(iParam+3),['\Delta' paramLabel{iParam}],'FontWeight','bold','FontSize',12);
    end
    if mod(iParam,3)==0 || mod(iParam,3)==1
        xlabel(hPlot(iParam),'Frequency (Hz)','FontWeight','bold','FontSize',12);
        xlabel(hPlot(iParam+3),'Frequency (Hz)','FontWeight','bold','FontSize',12);
    end

    yRange{iParam} = [min(yLim{iParam}) max(yLim{iParam})];
    yRange{iParam+3} = [min(yLim{iParam+3}) max(yLim{iParam+3})];

    makeShadedRegion(hPlot(iParam),freqRangeList{1},yRange{iParam},color1);
    makeShadedRegion(hPlot(iParam),freqRangeList{2},yRange{iParam},color2);

    makeShadedRegion(hPlot(iParam+3),freqRangeList{1},yRange{iParam+3},color1);
    makeShadedRegion(hPlot(iParam+3),freqRangeList{2},yRange{iParam+3},color2);

    if iParam==1
        legend(hPlot(iParam),'','','',[strList{1} '(' num2str(subjectNum{jDat}{1}) ')'],'',[strList{2} '(' num2str(subjectNum{jDat}{2}) ')'],'FontSize',10,'Location','northeast');
    end
end

%% slope vs HFD & power vs HFD
posSlope = ismember(centreFreqNew,freqValsHFD);
posHFD = ismember(freqValsHFD,centreFreqNew);

minFreqVal = 4; maxFreqVal = 950;
freqRangeWidth = 100;
freqValsPSD = freqValsHFD(posHFD);
for jFreq = 1:length(freqValsPSD)
    freq_range{jFreq} = [freqValsPSD(jFreq)-freqRangeWidth/2,freqValsPSD(jFreq)+freqRangeWidth/2];
    if freq_range{jFreq}(1)<minFreqVal
        freq_range{jFreq}(1) = minFreqVal;
    end

    if freq_range{jFreq}(2)>maxFreqVal
        freq_range{jFreq}(2) = maxFreqVal;
    end

    freqPos = freqVals{1}>=freq_range{jFreq}(1) & freqVals{1}<=freq_range{jFreq}(2);
    for iNum = 1:numGroups
        power{iNum}(:,jFreq) = log10(sum(10.^(dataPSD{iNum}(:,freqPos)),2));
    end
end

hPlot4 = subplot(numRows,2,7); box off;
hPlot5 = subplot(numRows,2,8); box off;
hPlot2 = subplot(numRows,2,9); box off;
hPlot3 = subplot(numRows,2,10); box off;

if medianFlag
    axes(hPlot2);
    for iNum = 1:numGroups
        plot(squeeze(nanmedian(dataParam{iNum}(:,posSlope),1)),squeeze(nanmedian(dataHFD{iNum}(:,posHFD),1)),'.','color',colorNames(iNum,:),'MarkerSize',10); hold on;
    end

    x1 = squeeze(nanmedian(dataParam{2}(:,posSlope),1)-nanmedian(dataParam{1}(:,posSlope),1));
    y1 = squeeze(nanmedian(dataHFD{2}(:,posHFD),1)-nanmedian(dataHFD{1}(:,posHFD),1));


    axes(hPlot4);
    for iNum = 1:numGroups
        plot(squeeze(nanmedian(power{iNum},1)),squeeze(nanmedian(dataHFD{iNum}(:,posHFD),1)),'.','color',colorNames(iNum,:),'MarkerSize',10); hold on;
    end

    x2 = squeeze(10*(nanmedian(power{2},1)-nanmedian(power{1},1)));
    y2 = squeeze(nanmedian(dataHFD{2}(:,posHFD),1)-nanmedian(dataHFD{1}(:,posHFD),1));

else
    for iNum = 1:numGroups
        plot(hPlot2,squeeze(nanmean(dataParam{iNum}(:,posSlope),2)),squeeze(nanmean(dataHFD{iNum}(:,posHFD),2)),'color',colorNames(iNum,:),'LineWidth',2); hold on;
    end

    axes(hPlot4);
    for iNum = 1:numGroups
        plot(squeeze(nanmean(power{iNum},1)),squeeze(nanmean(dataHFD{iNum}(:,posHFD),1)),'color',colorNames(iNum,:),'LineWidth',2); hold on;
    end

    x1 = squeeze(nanmean(dataParam{2}(:,posSlope),1)-nanmean(dataParam{1}(:,posSlope),1));
    y1 = squeeze(nanmean(dataHFD{2}(:,posHFD),1)-nanmean(dataHFD{1}(:,posHFD),1));

    x2 = squeeze(10*(nanmean(power{2},1)-nanmean(power{1},1)));
    y2 = squeeze(nanmean(dataHFD{2}(:,posHFD),1)-nanmean(dataHFD{1}(:,posHFD),1));

end
axes(hPlot3);
plot(x1,y1,'x','MarkerSize',6,'Color',displaySettingsDiff.colorNames); hold on;
mdl = fitlm(x1,y1);
yFit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*x1;
plot(x1,yFit,'LineWidth',2,'LineStyle','--','Color','k');
[C1,p1] = corrcoef(x1,y1);
text(-0.35,-0.002,['\rho = ' num2str(round(C1(1,2),2))]);
if p1<0.001
    text(-0.45,-0.004,['p = ' num2str(p1(1,2),'%.1e')]);
else
    text(-0.35,-0.004,['p = ' num2str(round(p1(1,2),3))]);
end

axes(hPlot5);
plot(x2,y2,'x','MarkerSize',6,'Color',displaySettingsDiff.colorNames); hold on;
mdl2 = fitlm(x2,y2);
yFit2 = mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2)*x2;
plot(x2,yFit2,'LineWidth',2,'LineStyle','--','Color','k');
[C2,p2] = corrcoef(x2,y2);
text(-0.45,-0.002,['\rho = ' num2str(round(C2(1,2),2))]);
if p2<0.001
    text(-0.45,-0.004,['p = ' num2str(p2(1,2),'%.1e')]);
else
    text(-0.45,-0.004,['p = ' num2str(round(p2(1,2),3))]);
end

set(hPlot2,'FontWeight','Bold','LineWidth',1,'FontSize',10,'box','off');
set(hPlot3,'FontWeight','Bold','LineWidth',1,'FontSize',10,'box','off');
set(hPlot4,'FontWeight','Bold','LineWidth',1,'FontSize',10,'box','off');
set(hPlot5,'FontWeight','Bold','LineWidth',1,'FontSize',10,'box','off');

xlabel(hPlot2,'Slope','FontSize',12,'FontWeight','bold');
ylabel(hPlot2,'HFD','FontSize',12,'FontWeight','bold');
xlabel(hPlot3,'\Delta Slope','FontSize',12,'FontWeight','bold');
ylabel(hPlot3,'\Delta HFD','FontSize',12,'FontWeight','bold');

xlabel(hPlot4,'log_{10}(Power(\muV^2))','FontSize',12,'FontWeight','bold');
ylabel(hPlot4,'HFD','FontSize',12,'FontWeight','bold');
xlabel(hPlot5,'\Delta Power (dB)','FontSize',12,'FontWeight','bold');
ylabel(hPlot5,'\Delta HFD','FontSize',12,'FontWeight','bold');


function makeShadedRegion(hPlot,x,y,color,faceAlphaVal)
if ~exist('faceAlphaVal','var');  faceAlphaVal=0.3;   end
xx = [x(1) x(1) x(2) x(2)];
yy = [y(1)  y(2) y(2) y(1)];
patch(hPlot,xx,yy,color,'LineStyle','none','FaceAlpha',faceAlphaVal);
h = get(hPlot,'Children');
set(hPlot,'Children',circshift(h,length(h)-1));
end
