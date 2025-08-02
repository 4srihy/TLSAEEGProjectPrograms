%% it generates fig 6 contain power and hfd matched for all frequency ranges
plotFigures;
%plotting TF Power and HFD
clear powerBl powerSt tfPower dim0 dimS0 dimDiff dimTF
bonferroniFlag = 0;
titleMatch = [{'Alpha'} {'Slow Gamma'} {'Fast Gamma'}];
%1:alpha, 2: slow gamma , 3: fastGamma
ElecGroupToPlot = 6;
medianFlag=1;
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
            % tfPower{i}(j,:,:) = x.tfPower;
            powerBl{i}(j,:) = x.blPowerVsFreq;
            powerSt{i}(j,:) = x.stPowerVsFreq;
        end
    end
    deltaPower0{i} = 10*(log10(powerSt{i})-log10(powerBl{i}));
end

freqValsLinearPower = x.freqVals;

%% Loading HFD data
folderName = fullfile(folderSourceString{1},projectName,'higFracDim');
stFlag = 1;

clear dim0 HFD HFDBl HFDSt diffHFDSt diffHFDBL diffHFD
for i=1:length(subjectNameList)
    numHFD(i) = 0;% for linear hfd file
    if diseaseFlag
        if i==1;     medianFlag=1;    else;    medianFlag = 0; end
    end
    for j=1:length(subjectNameList{i})
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
                    if stFlag;  dimS0{i}(j,k,:) = squeeze(nanmean(xj1.fracDim.higSt2{k}(:,1,:),3));   end
                end
            end
        end
    end
    dim0{i}(dim0{i}==0) = nan;
    if stFlag
        dimS0{i}(dimS0{i}==0) = nan;
        dimDiff{i} = dimS0{i}-dim0{i};
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

for iNum=1: numGroup
    dataHFD0{iNum} = squeeze(nanmean(dimDiff{iNum}(:,:,electrodeGroupList{ElecGroupToPlot}),3)); %for linear plot
end

capType = 'actiCap64';
kPos=[];

%% Matching

totIter=10;
powerMatchingBinWidth = 3;
matchedSubjectNameList = getPowerMatchedSubjects(subjectNameList,folderName1,'DiffdB',totIter,powerMatchingBinWidth);
figure;
numRows = 2;
numCols = 3;
freqRanges = {[8 12], [20 36], [40 66]};

for numfreqRangeToMatch = 1:3
    clear kPos numHFD numPower
    for iTer = 1:totIter
        for iNum = 1:length(matchedSubjectNameList{numfreqRangeToMatch}{iTer})
            for iMatch = 1:length(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum})
                for iSub = 1:length(subjectNameList{iNum})
                    if strcmp(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum}{iMatch},subjectNameList{iNum}{iSub})
                        kPos{iNum}{iTer}(iMatch) = iSub;
                        continue;
                    end
                end

            end
            numHFD(iNum) = length(matchedSubjectNameList{numfreqRangeToMatch}{iTer}{iNum});
            numPower(iNum) = numHFD(iNum);
        end
    end


    for i=1:6;       jh{i} = ':';    end % for size balancing

    clear HFD dtfPowerMean dataHFD1 deltaPower1
    for iNum=1:numGroup
        for iTer = 1:totIter
            dataHFD1{iNum}(iTer,jh{:}) = dataHFD0{iNum}(kPos{iNum}{iTer},jh{:});
            deltaPower1{iNum}(iTer,jh{:}) = deltaPower0{iNum}(kPos{iNum}{iTer},jh{:});
        end

        dataHFD{iNum} = squeeze(nanmean(dataHFD1{iNum}(jh{:}),1));
        deltaPower{iNum} = squeeze(nanmean(deltaPower1{iNum}(jh{:}),1));
    end


    clear parameter

    parameter{7} = deltaPower;
    parameter{8} = dataHFD;

    labelList = [{'Mid'} {'Old'} {'Diff'}];


    displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
    xScaleVar = 'linear';
    colorNames = [1 0 0; 0 0 1; 1 0 1;0.49 0.18 0.56; 0.5 0.5 0.5; 0 1 1 ; 0 1 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%hot(8); %colorNames([1:3,end-2:end],:) = [];
    displaySettings.colorNames = colorNames;
    displaySignificanceFlag=1;

    param{1} = parameter{7};    xVals{1} = freqValsLinearPower;     paramLabel{1} = 'Power (dB)';   numParam{1}  = numPower;   yLim{1} = [-1 1.5];
    param{2} = parameter{8};    xVals{2} = freqValsHFDLinear;       paramLabel{2} = 'HFD';          numParam{2} = numHFD;       yLim{2} = [-2*10^-4 3*10^-4];

    for iParam = 1:2
        hplotLinear = subplot(numRows,numCols,(iParam-1)*3+numfreqRangeToMatch);
        displayAndcompareData(hplotLinear,param{iParam},xVals{iParam},displaySettings,yLim{iParam},displaySignificanceFlag,medianFlag,'','','',bonferroniFlag); xlim([6 100]);
        hold on;
        if numfreqRangeToMatch==1
            ylabel(hplotLinear,['\Delta' paramLabel{iParam}],'fontsize',12,'FontWeight','bold');
        end
        if iParam==2;    xlabel('Frequency (Hz)','fontsize',12,'FontWeight','bold'); end

        makeBox(hplotLinear,freqRanges{numfreqRangeToMatch},hplotLinear.YLim,'k',1,'-','V');

        plot(hplotLinear,hplotLinear.XLim,0.*hplotLinear.XLim,'k','LineWidth',0.8);
        set(hplotLinear,'LineWidth',1,'FontSize',12,'FontWeight','bold');

        if iParam==1
            title([titleMatch{numfreqRangeToMatch} ' Power Matched'],'fontsize',14,'FontWeight','bold')
            [l1,l2]=legend(hplotLinear,'',[strList{1} '(' num2str(numParam{iParam}(1)) ')'],'',[strList{2} '(' num2str(numParam{iParam}(2)) ')'],'Box','off');

            l2(3).XData(1) = 0.15;
            l2(5).XData(1) = 0.15;
        end
    end
end

a = annotation('textbox',[.058 0.9415 0.0228 0.0579],'String','(A)','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.3509 0.9415 0.0228 0.0579],'String','(B)','FontSize',16,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.6322 0.9415 0.0228 0.0579],'String','(C)','FontSize',16,'EdgeColor','none','FontWeight','bold');
