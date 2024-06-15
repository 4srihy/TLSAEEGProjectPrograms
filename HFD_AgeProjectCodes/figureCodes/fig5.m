%plotting TF Power and HFD
clear powerBl powerSt tfPower dim0 dimS0 dimDiff dimTF 
bonferroniFlag = 0;
matchingFlag = 0;
numfreqRangeToMatch = 1; %1:alpha, 2: slow gamma , 3: fastGamma
ElecGroupToPlot = 6;
medianFlag=0;
freqWidthHfd = 12;
diseaseFlag = 0;

%% Loading power data
folderName1 =  fullfile(folderSourceString{1},projectName,'SF_ORI');
blRange = [-0.5 0];
timeLims = [-0.25 1.02];
for i=1:length(subjectNameList) 
    numPower(i)=0;
    for j=1:length(subjectNameList{i})
        subjectName = subjectNameList{i}{j};
        fileName = fullfile(folderName1,[subjectName '_' refType '_stRange_250_750.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectName ' does not exist']);
        else
            numPower(i) = numPower(i)+1;
            x = load(fileName);
            tfPower{i}(j,:,:) = x.tfPower;
            powerBl{i}(j,:) = x.blPowerVsFreq;
            powerSt{i}(j,:) = x.stPowerVsFreq;
        end
    end
    deltaPower{i} = 10*(log10(powerSt{i})-log10(powerBl{i}));
    timesTFPos = x.timeValsTF>=timeLims(1) & x.timeValsTF<=timeLims(2);
    timeValsTFPower = x.timeValsTF(timesTFPos);
    blPosTF = x.timeValsTF>=blRange(1) & x.timeValsTF<=blRange(2);
    TFPowerBL{i} = nanmean(tfPower{i}(:,blPosTF,:),2);
    dtfPower{i} = 10.*log10(tfPower{i}(:,timesTFPos,:)./TFPowerBL{i});
    
end
freqValsTFPower = x.freqValsTF;
freqValsLinearPower = x.freqVals;

%% Loading HFD data
folderName = fullfile(folderSourceString{1},projectName,'higFracDim');
stFlag = 1;

clear dim0 HFD HFDBl HFDSt diffHFDSt diffHFDBL diffHFD
for i=1:length(subjectNameList)
    numHFD(i) = 0;% for linear hfd file
    numHFDTF(i) = 0; %for hfd tf file
    if diseaseFlag
        if i==1;     medianFlag=1;    else;    medianFlag = 0; end
    end
    for j=1:length(subjectNameList{i})
        fileName = fullfile(folderName,[subjectNameList{i}{j} '_' refType '_TF_nonoverlapping_kmax5.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
        else
            numHFDTF(i) = numHFDTF(i)+1;
            xj = load(fileName);
            % we want to create a mtrix of (subject X elec X time X freq)
            for ktime=1:length(xj.fracDim.higBL2)
                for kFreq = 1:length(xj.fracDim.higBL2{ktime})
                    z1 = cell2mat(xj.fracDim.higBL2{ktime}(1,kFreq));
                    for ci = 1:size(xj.fracDim.higBL2{ktime},1)-1
                        z1 = cat(3,z1,cell2mat(xj.fracDim.higBL2{ktime}(ci+1,kFreq)));
                    end
                    dim0TF{i}(j,:,ktime,kFreq) = squeeze(nanmean(z1(:,1,:),3));
                end
            end
        end

        fileName2 = fullfile(folderName,[subjectNameList{i}{j} '_' refType '_freqWidth' num2str(freqWidthHfd) '_kmax5.mat']);
        if ~isfile(fileName2)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
        else
            numHFD(i) = numHFD(i)+1;
            xj1 = load(fileName2);
            for k=1:length(xj1.fracDim.higBL2)
                if freqWidthHfd==24
                    zc  = bipolarMean(xj1.fracDim.higBL2{k}',0);
                    dim0{i}(j,k,:) = zc(1,:);
                    if stFlag
                        zcs = bipolarMean(xj1.fracDim.higSt2{k}',0);
                        dimS0{i}(j,k,:) = zcs(1,:);
                    end
                else
                    dim0{i}(j,k,:) = squeeze(nanmean(xj1.fracDim.higBL2{k}(:,1,:),3));
                    if stFlag  dimS0{i}(j,k,:) = squeeze(nanmean(xj1.fracDim.higSt2{k}(:,1,:),3));   end
                end
            end
        end
    end
    dim0{i}(dim0{i}==0) = nan;
    if stFlag
        dimS0{i}(dimS0{i}==0) = nan;
        dimDiff{i} = dimS0{i}-dim0{i};
    end

    if numHFD(i)~=numHFDTF(i)
        disp('no. of subjects for linear and tf HFD plots are not same');
    end
end

%% electrode group
if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
    [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    groupNameList{length(groupNameList)+1} = 'highPriority';
end

numGroup = length(subjectNameList);
%linear HFD
if freqWidthHfd == 24
    freqValsHFDLinear = xj1.fracDim.centreFreq(1:29);
else
    freqValsHFDLinear = squeeze(xj1.fracDim.centreFreq(:,:,1));
end


%TF HFD
timeVar  =xj.fracDim.timeVar(1,1:end-1);
freqVar = xj.fracDim.freqVar(1,1:end-1);
for iNum=1: numGroup
     HFD0{iNum} = squeeze(nanmean(dim0TF{iNum}(:,electrodeGroupList{ElecGroupToPlot},:,:),2));
     dataHFD0{iNum} = squeeze(nanmean(dimDiff{iNum}(:,:,electrodeGroupList{ElecGroupToPlot}),3)); %for linear plot
end

numBl = find(timeVar<0,1,'last');

capType = 'actiCap64';
kPos=[];

%% Matching
if matchingFlag
    totIter=10;
    powerMatchingBinWidth = 3;
    matchedSubjectNameList = getPowerMatchedSubjects(subjectNameList,folderName1,'DiffdB',totIter,powerMatchingBinWidth);
    for iTer = 1:totIter
        for iNum = 1:length(matchedSubjectNameList{numfreqRangeToMatch}{iTer})
            for iMatch = 1:length(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum})
                for iSub = 1:length(subjectNameList{iNum})
                    if strcmp(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum}{iMatch},subjectNameList{iNum}{iSub})
                        kPos{iNum}{iTer}(iMatch) = iSub;
                        continue;
                    end
                end
                numHFD(iNum) = length(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum});
                numPower(iNum) = numHFD(iNum);
            end

        end
    end
else
    totIter=1;
    kPos{1}{1} = ':';   kPos{2}{1} = ':';
end

for i=1:6;       jh{i} = ':';    end % for size balancing

clear HFD dtfPowerMean dataHFD1 deltaPower1
for iNum=1:numGroup
   for iTer = 1:totIter
    if medianFlag
        HFD{iNum}(iTer,jh{:}) = squeeze(nanmedian(HFD0{iNum}(kPos{iNum}{iTer},jh{:}),1));
        dtfPowerMean{iNum}(iTer,jh{:}) = squeeze(nanmedian(dtfPower{iNum}(kPos{iNum}{iTer},jh{:}),1));
    else
        HFD{iNum}(iTer,jh{:}) = squeeze(nanmean(HFD0{iNum}(kPos{iNum}{iTer},jh{:}),1));
        dtfPowerMean{iNum}(iTer,jh{:}) = squeeze(nanmean(dtfPower{iNum}(kPos{iNum}{iTer},jh{:}),1));
    end
    dataHFD1{iNum}(iTer,jh{:}) = dataHFD0{iNum}(kPos{iNum}{iTer},jh{:});
    deltaPower1{iNum}(iTer,jh{:}) = deltaPower{iNum}(kPos{iNum}{iTer},jh{:});
   end

    HFD{iNum} = squeeze(nanmean(HFD{iNum}(jh{:}),1));
    dtfPowerMean{iNum} = squeeze(nanmean(dtfPowerMean{iNum}(jh{:}),1));
    HFDBl{iNum} = nanmean(HFD{iNum}(1:numBl,jh{:}),1);
    HFDSt{iNum} = HFD{iNum}-HFDBl{iNum};
    dataHFD{iNum} = squeeze(nanmean(dataHFD1{iNum}(jh{:}),1));
   deltaPower{iNum} = squeeze(nanmean(deltaPower1{iNum}(jh{:}),1));
end


clear parameter
parameter{1} = dtfPowerMean{1};
parameter{2} = dtfPowerMean{2};
parameter{3} = dtfPowerMean{2}-dtfPowerMean{1};

parameter{4} = HFDSt{1};
parameter{5} = HFDSt{2};
parameter{6} = HFDSt{2}-HFDSt{1};

parameter{7} = deltaPower;
parameter{8} = dataHFD;

labelList = [{'Mid'} {'Old'} {'Old-Mid'}];


%% plotting
figure;
numRows = 2;
numCols = 4;
freqRanges = {[8 12], [20 36], [40 66]};
for i=1:6
    hplot = subplot(numRows,numCols,floor((i-0.1)/3)+i+1);
    if i<4
        timePlot = timeValsTFPower;
        pcolor(hplot,timeValsTFPower,freqValsTFPower,parameter{i}');
        title(labelList{i},'fontsize',13);
        colormap(hplot,'jet');
    else
        timePlot = timeVar;
        pcolor(hplot,timeVar,freqVar,parameter{i}(1:length(timeVar),:)');
        d1 = colormap(hplot,'jet');
        d2 = flipud(d1);
        colormap(hplot,d2);
        xlabel(hplot,'Time (s)','fontsize',14,'FontWeight','bold');
    end
    shading('interp')
    %caxis([-2*10^-4,2*10^-4])

    if i<3
        caxis([-1.5 1.5]);
    elseif i==3
        caxis([-0.5 0.5]);
    elseif i<6
        caxis([-2*10^-4,2*10^-4])
    else
        caxis([-0.8*10^-4,0.8*10^-4])
    end
    c = colorbar;

    if i==3
        c.Label.String = '\DeltaPower (dB)';
        c.Label.FontSize = 13;
        c.Label.FontWeight = 'bold';
    elseif i==6
        c.Label.String = '\DeltaHFD';
        c.Label.FontSize = 13;
        c.Label.FontWeight = 'bold';
    end

    ylim([0 100]);
    makeBox(hplot,[timePlot(1) timePlot(end)],freqRanges{1},'w',2,'-','H');
    makeBox(hplot,[timePlot(1) timePlot(end)],freqRanges{2},'w',2,'--','H');
    makeBox(hplot,[timePlot(1) timePlot(end)],freqRanges{3},'w',2,':','H');
   
    if ismember(i,[1 4])
        ylabel(hplot,'Frequency (Hz)','fontsize',14,'FontWeight','bold');
    end
    set(hplot,'FontSize',12,'FontWeight','Bold');
end

displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
xScaleVar = 'linear';
colorNames = [1 0 0; 0 0 1; 1 0 1;0.49 0.18 0.56; 0.5 0.5 0.5; 0 1 1 ; 0 1 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%hot(8); %colorNames([1:3,end-2:end],:) = [];
displaySettings.colorNames = colorNames;
displaySignificanceFlag=1;

param{1} = parameter{7};         xVals{1} = freqValsLinearPower;     paramLabel{1} = 'Power (dB)';   numParam{1}  = numPower;   yLim{1} = [-1 1.5];
param{2} = parameter{8};         xVals{2} = freqValsHFDLinear;       paramLabel{2} = 'HFD';          numParam{2} = numHFD;      yLim{2} = [-2*10^-4 3*10^-4];    

for iParam = 1:2
    %hplotLinear = subplot(numRows,numCols,numCols*(iParam-1)+1);
    hplotLinear = subplot('Position',[0.0925 0.5840-(iParam-1)*0.4644 0.1718 0.3513]);
    displayAndcompareData(hplotLinear,param{iParam},xVals{iParam},displaySettings,yLim{iParam},displaySignificanceFlag,medianFlag,'','','',bonferroniFlag); xlim([6 100]);
     hold on;
    ylabel(hplotLinear,['\Delta' paramLabel{iParam}]);
    if iParam==2;    xlabel('Frequency (Hz)'); end
    
    if matchingFlag
    makeBox(hplotLinear,freqRanges{numfreqRangeToMatch},yLim{iParam},'k',1,'-','V');
    else
       for iFreq = 1:length(freqRanges)
            makeBox(hplotLinear,freqRanges{iFreq},yLim{iParam},colorNames(iFreq+2,:),1,'-','V');
       end
    end
   plot(hplotLinear,hplotLinear.XLim,0.*hplotLinear.XLim,'k','LineWidth',0.8);
   set(hplotLinear,'LineWidth',1,'FontSize',12,'FontWeight','Bold');
   
    if iParam==1
        [l1,l2] = legend(hplotLinear,'',strList{1},'',strList{2},'Box','off','FontSize',12);
        set(l2,'lineWidth',2)
        l2(3).XData(1) = 0.4;
        l2(5).XData(1) = 0.4;
        %legend(hplotLinear,'',[strList{1} '(' num2str(numParam{iParam}(1)) ')'],'',[strList{2} '(' num2str(numParam{iParam}(2)) ')'],'Box','off','FontSize',12);
    end

end

a = annotation('textbox',[.0158 0.9415 0.0228 0.0579],'String','(A)','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.2713 0.9415 0.0228 0.0579],'String','(B)','FontSize',16,'EdgeColor','none','FontWeight','bold');
