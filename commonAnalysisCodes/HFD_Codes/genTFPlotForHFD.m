function fracDim = genTFPlotForHFD(eegData,timeVals,electrodeList,refType,timeBin,freqBin,maxFreq, plotFlag, plotType, hPlot,minTime, maxTime)
if ~exist('plotType','var');      plotType = 'diff';   end %could be 'diff','abs','both'
if ~exist('refType','var');      refType = 'unipolar';   end
if ~exist('timeBin','var');      timeBin = 100;   end
if ~exist('freqBin','var');      freqBin = 8;   end
if ~exist('maxFreq','var');     maxFreq = 151;   end
if ~exist("plotFlag","var");     plotFlag = 0;   end
if ~exist('minTime','var');     minTime = -300;     end %in ms for TLSA
if ~exist("maxTime",'var');     maxTime = 1200;     end  % for TLSA

%if ~exist('kMax','var')         kMax  = 20;         end

minFreqVal = 1;
% minTime = -300;    %in ms
% maxTime = 1200;

freqVar = minFreqVal:freqBin:maxFreq;
timeVar = (minTime+timeBin/2:timeBin:maxTime-timeBin/2)*10^-3;
minTime = minTime*10^-3;    %in ms
maxTime = maxTime*10^-3;

% for overlapping time windows
% timeRangeWidth = timeRangeWidth*10^-3;
%
% for jTime = 1:length(timeVar)
%     timeRanges{jTime} = [timeVar(jTime)-timeRangeWidth/2,timeVar(jTime)+timeRangeWidth/2];
%     if timeRanges{jTime}(1)<minTime
%         timeRanges{jTime}(1) = minTime;
%     end
%
%     if timeRanges{jTime}(2)>maxTime
%         timeRanges{jTime}(2) = maxTime;
%     end
% end

%for non overlapping
for jTime = 1:length(timeVar)-1
    timeRanges{jTime}(1) = timeVar(jTime);
    timeRanges{jTime}(2) = timeVar(jTime+1);
end

for jFreq = 1:length(freqVar)-1
    bandPassFilterFreq{jFreq}(1) = freqVar(jFreq);
    bandPassFilterFreq{jFreq}(2) = freqVar(jFreq+1);
end

%
% maxFreqVal = 200;
% centreFreq =12:10:192;%12:10:142;
% freqRangeWidth = 24;
% clear bandPassFilterFreq freqLength

awayFromCentre = 1;
for freqAvoidcentre = 50:50:250
    stopFreqRange{freqAvoidcentre/50} = [freqAvoidcentre-awayFromCentre freqAvoidcentre+awayFromCentre];
end

%bandPassFilterFreq = {[1 90]};% {[8 12], [20 35],[60 90], [110 140], [160 190],[210 240],[260 290]};
%bandPassFilterFreq = {[7 13], [20 35],[60 90]}; % alpha, slow gamma, and aperiodic
%timeRanges = {[-0.5 -0.25],[-0.25 0], [0 .25], [0.25 0.5],[0.5 0.75]};
%timeRanges = {[-0.5 -0.4],[-0.5 -0.3], [-0.5 -0.2], [-0.5 -0.1],[-0.5 -0.01]};
%timeRanges = {[-0.5 -0.25],[0.5 0.75]};
numTimeRanges = length(timeRanges);
numFreqRanges = length(bandPassFilterFreq);

rmseMax = 0.05;
fastFlag = 1;
optimiseRmseFlag = 0;
%kMax =[15	15	15	15	15	15	12	10	9	8	8	7	6	6	5	5	4	4	4	4	4];%15*ones(1,numFreqRanges);%10*ones(1,5);%[10 11 12 13 14 15];%[10 10 10 7 5 4 3]; % for 1-90 :kmax=10
kMax = 5*ones(1,numFreqRanges);

% Get good positions
Fs = round(1/(timeVals(2)-timeVals(1)));

for iRange = 1:numTimeRanges
    range = timeRanges{iRange};
    rangePos = round(diff(range)*Fs);
    timePos{iRange} = find(timeVals>=range(1),1)+ (1:rangePos);
end

% Initialize
numElectrodes = length(electrodeList);
numSides = 1;%length(electrodeList);
numTrials = size(eegData,2);

numAnalysedElecs = repmat(numElectrodes,numSides,1);
hW = waitbar(0,'Analysing electrodes...');

%for average refernce scheme
if strcmp(refType,'average')
    allData = [];
    for iElec=1:numElectrodes
        analogData = squeeze(eegData(iElec,:,:));
        allData = cat(3,allData,analogData);
    end
    meanElecData = squeeze(nanmean(allData,3));
end

%filters
d{1} = designfilt('highpassiir','FilterOrder',8, 'PassbandFrequency',1,'PassbandRipple',0.2,'SampleRate',Fs);
for iFilt = 1:length(bandPassFilterFreq)
    b{iFilt} = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',bandPassFilterFreq{iFilt}(1),'PassbandFrequency2',bandPassFilterFreq{iFilt}(2),'PassbandRipple',0.2,'SampleRate',Fs);
end
tic;
for iRange = 1:numTimeRanges
    numTimes = length(timePos{iRange});
    iEs = 0;

    for iElec=1:numElectrodes % For each electrode or electrode pair

        for iSide = 1:numSides

            if ~isnan(eegData(iElec,:,:))
                iEs = iEs+1;
                hW = waitbar(iEs/(numElectrodes*numSides),hW,'Analysing electrodes...');

                if strcmp(refType,'unipolar') %length(electrodeList{iSide}{iElec})==1 % Single Electrode
                    eeg = squeeze(eegData(electrodeList(iElec),:,:));
                elseif strcmp(refType,'average')
                    eeg = squeeze(eegData(electrodeList(iElec),:,:)) - meanElecData;
                elseif strcmp(refType,'bipolar')% for bipolar referencing
                    chan1 = electrodeList{iSide}{iElec}(1);
                    chan2 = electrodeList{iSide}{iElec}(2);
                    eeg = squeeze(eegData(chan1,:,:) - eegData(chan2,:,:));
                end

                timeData{iRange}(iSide,iElec,:,:) =  eeg(:,timePos{iRange});

                %Filtered data
                dataNewBl{iRange}(iSide,iElec,:,:) = filtfilt(d{1},squeeze(timeData{iRange}(iSide,iElec,:,:))');
                for iFilt = 1:length(bandPassFilterFreq)
                    dataNewBandBl{iRange}{iFilt}(iSide,iElec,:,:)   = filtfilt(b{iFilt},squeeze(timeData{iRange}(iSide,iElec,:,:))');
                end

                for iTrial = 1:numTrials
                    for iFilt = 1:length(bandPassFilterFreq)
                        [higBL2{iRange}{iFilt}(iSide,iElec,iTrial,1),higBL2{iRange}{iFilt}(iSide,iElec,iTrial,2),higBL2{iRange}{iFilt}(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBandBl{iRange}{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);
                        if optimiseRmseFlag
                            kMax(iFilt) = higBL2{iRange}{iFilt}(iSide,iElec,iTrial,3);
                            kMax(iFilt+1) = kMax(iFilt);
                        end
                    end
                    optimiseRmseFlag = 0;
                    fastFlag = 1;
                end

            else

                timeData{iRange}(iSide,iElec,[1:numTrials],[1:numTimes]) = nan;
                dataNewBl{iRange}(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                for iFilt = 1:length(bandPassFilterFreq)
                    higBL2{iRange}{iFilt}(iSide,iElec,[1:numTrials],[1:3]) = nan;
                    dataNewBandBl{iRange}{iFilt}(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                end

            end
        end

    end

    %   close(hW);

    %% Taking mean across trials
    dataNewBl{iRange} = (squeeze(nanmean(dataNewBl{iRange},4)))';

    for iFilt = 1:length(bandPassFilterFreq)
        dataNewBandBl{iRange}{iFilt} = (squeeze(nanmean(dataNewBandBl{iRange}{iFilt},4)))';
        higBL2{iRange}{iFilt} = squeeze(nanmean(higBL2{iRange}{iFilt},3));
    end

end
toc;
%% HFD for TF
if plotFlag

    for iRange = 1:numTimeRanges
        for iFreq = 1:numFreqRanges
            HFD(iRange,iFreq) = nanmean(higBL2{iRange}{iFreq}(1));

        end
    end
    %% diff HFD
    numBl = max(find(timeVar<0));
    HFDBl = mean(HFD(1:numBl,:),1);
    diffHFD = HFD-HFDBl;

    %% Plot

    switch plotType
        case 'abs'
            axes(hPlot);
            pcolor(timeVar(1:numTimeRanges),freqVar(1:numFreqRanges),HFD')
            colorbar;   shading('interp');
            xlabel('time(s)');
            ylabel('frequency (Hz)')
            title('HFD')
        case 'diff'
            pcolor(hPlot,timeVar(1:numTimeRanges),freqVar(1:numFreqRanges),diffHFD')
            colorbar;   shading('interp');
            caxis([-0.0025 0.0025]);
            title('\Delta HFD')

        case 'both'
            pcolor(hPlot(1),timeVar(1:numTimeRanges),freqVar(1:numFreqRanges),HFD')
            colorbar(hPlot(1));   shading(hPlot(1),'interp');
            xlabel('time(s)');
            ylabel('frequency (Hz)')
            title('HFD');
            pcolor(hPlot(2),timeVar(1:numTimeRanges),freqVar(1:numFreqRanges),diffHFD')
            colorbar(hPlot(2));   shading(hPlot(2),'interp');
            caxis(hPlot(2),[-0.0005 0.0005]);
            title('\Delta HFD')
    end
    colormap('jet')
    xlabel('Time(s)');
    ylabel('Frequency (Hz)')

    fracDim.HFD  = HFD;
    fracDim.diffHFD = diffHFD;
end

fracDim.kmax = kMax;
fracDim.higBL2 = higBL2;
fracDim.bandPassFilterFreq = bandPassFilterFreq;

fracDim.timeRanges = timeRanges;
fracDim.timeVar = timeVar;
fracDim.freqVar = freqVar;

end