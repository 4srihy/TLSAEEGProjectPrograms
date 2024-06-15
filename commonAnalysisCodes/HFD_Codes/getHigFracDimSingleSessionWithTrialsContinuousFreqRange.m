
function [data,fracDim,numAnalysedElecs]=getHigFracDimSingleSessionWithTrialsContinuousFreqRange(eegData,timeVals,electrodeList,stRange,stFlag,refType)

%minFreqVal = 1; centreFreq = 7:12:235; maxFreqVal  =250; freqRangeWidth = 12;
%minFreqVal = 4;
%centreFreq =12:10:82;  maxFreqVal = 90;
%centreFreq = 12:10:292; freqRangeWidth = 24;
%centreFreq = 32:10:292;  freqRangeWidth = 52; maxFreqVal = 340;

%for TLSA project
%centreFreq = 40:10:300;  freqRangeWidth = 100;maxFreqVal = 360;

%for meditation project
%centreFreq = 40:10:180;  freqRangeWidth = 100;maxFreqVal = 220;
minFreqVal = 1; centreFreq = 7:12:199; maxFreqVal  =210; freqRangeWidth = 12;

clear bandPassFilterFreq freqLength
for jFreq = 1:length(centreFreq)
    bandPassFilterFreq{jFreq} = [centreFreq(jFreq)-freqRangeWidth/2,centreFreq(jFreq)+freqRangeWidth/2];
    if bandPassFilterFreq{jFreq}(1)<minFreqVal
        bandPassFilterFreq{jFreq}(1) = minFreqVal;
    end

    if bandPassFilterFreq{jFreq}(2)>maxFreqVal
        bandPassFilterFreq{jFreq}(2) = maxFreqVal;
    end
end

awayFromCentre = 2;
for freqAvoidcentre = 50:50:250
    stopFreqRange{freqAvoidcentre/50} = [freqAvoidcentre-awayFromCentre freqAvoidcentre+awayFromCentre];
end
%if ~exist('kMax','var')         kMax  = 20;         end
%bandPassFilterFreq = {[1 90]};% {[8 12], [20 35],[60 90], [110 140], [160 190],[210 240],[260 290]};
%bandPassFilterFreq = {[7 13], [20 35],[60 90]}; % alpha, slow gamma, and aperiodic
if stFlag
    blRange = [-diff(stRange) 0];
    %kMax = 10; % for low resolution data
    fastFlag = 1; %no rmse data for fitting
else
    blRange = [timeVals(1)  timeVals(length(timeVals))];     stRange = blRange;
    fastFlag = 0; % rmse data for fitting exists and it is slow
    %kMax = 10;
end
rmseMax = 0.05;
fastFlag = 1;
optimiseRmseFlag = 0;
%kMax = [10	10	10	10	10	10	9	9	8	7	7	6	6	5	5	5	5	4	4	4	4	4	3	3	3	3	3	3]; %for freqWidth 52
%kMax = [10	10	10	10	10	10	10	10	10	9	8	8	7	6	6	6	5	5	5	4	4	4	4	3	3	3	3	3	3 3]; for freqwidth 24
kMax = 5*ones(1,length(bandPassFilterFreq));%[10 10 10 7 5 4 3]; % for 1-90 :kmax=10
% Get good positions
Fs = round(1/(timeVals(2)-timeVals(1)));
if round(diff(blRange)*Fs) ~= round(diff(stRange)*Fs)
    disp('baseline and stimulus ranges are not the same');
else
    range = blRange;
    rangePos = round(diff(range)*Fs);
    blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
    stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);
end

% Initialize
numElectrodes = length(electrodeList{1});
numSides = length(electrodeList);
numTrials = size(eegData,2);
numTimes = length(blPos);
numAnalysedElecs = repmat(numElectrodes,numSides,1);
stPowerVsFreq = [];
% hW = waitbar(0,'Analysing electrodes...'); iEs = 0;

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

    clear bStop
    iStop0 = 0;
    for iStop = 1:length(stopFreqRange)
        if (stopFreqRange{iStop}(1)<bandPassFilterFreq{iFilt}(2)  && stopFreqRange{iStop}(1)>bandPassFilterFreq{iFilt}(1)) || (stopFreqRange{iStop}(2)<bandPassFilterFreq{iFilt}(2)  && stopFreqRange{iStop}(2)>bandPassFilterFreq{iFilt}(1))
            iStop0 = iStop0+1;
            bStop{iStop0} = designfilt('bandstopiir','FilterOrder',4,'PassbandFrequency1',stopFreqRange{iStop}(1),'PassbandFrequency2',stopFreqRange{iStop}(2),'PassbandRipple',0.2,'SampleRate',2500);
        end
    end

    for iElec=1:numElectrodes % For each electrode or electrode pair
        %tic;
        for iSide = 1:numSides
            if ~isnan(eegData(electrodeList{iSide}{iElec},:,:))

                %         iEs = iEs+1;
                %         hW = waitbar(iEs/(numElectrodes*numSides),hW,'Analysing electrodes...');

                if strcmp(refType,'unipolar') %length(electrodeList{iSide}{iElec})==1 % Single Electrode
                    eeg = squeeze(eegData(electrodeList{iSide}{iElec},:,:));
                elseif strcmp(refType,'average')
                    eeg = squeeze(eegData(electrodeList{iSide}{iElec},:,:)) - meanElecData;
                elseif strcmp(refType,'bipolar')% for bipolar referencing
                    chan1 = electrodeList{iSide}{iElec}(1);
                    chan2 = electrodeList{iSide}{iElec}(2);
                    eeg = squeeze(eegData(chan1,:,:) - eegData(chan2,:,:));
                end

                blData(iSide,iElec,:,:) =  eeg(:,blPos);

                %Filtered data
                % dataNewBl(iSide,iElec,:,:) = filtfilt(d{1},squeeze(blData(iSide,iElec,:,:))');
                %  for iFilt = 1:length(bandPassFilterFreq)
                dataNewBandBl{iFilt}(iSide,iElec,:,:)   = filtfilt(b{iFilt},squeeze(blData(iSide,iElec,:,:))');
                % end

                if exist('bStop','var')
                    for iStopFilt = 1:length(bStop)
                        dataNewBandBl{iFilt}(iSide,iElec,:,:) = filtfilt(bStop{iStopFilt},squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,:)));
                    end
                end

                if stFlag
                    stData(iSide,iElec,:,:) =  eeg(:,stPos);
                    dataNewSt(iSide,iElec,:,:) = filtfilt(d{1},squeeze(stData(iSide,iElec,:,:))');
                    % for iFilt = 1:length(bandPassFilterFreq)
                    dataNewBandSt{iFilt}(iSide,iElec,:,:)   = filtfilt(b{iFilt},squeeze(stData(iSide,iElec,:,:))');
                    %end
                    if exist('bStop','var')
                        for iStopFilt = 1:length(bStop)
                            dataNewBandSt{iFilt}(iSide,iElec,:,:) = filtfilt(bStop{iStopFilt},squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,:)));
                        end
                    end
                end

                for iTrial = 1:numTrials
                    %[higBLOr(iSide,iElec,iTrial,1),higBLOr(iSide,iElec,iTrial,2),higBLOr(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(blData(iSide,iElec,:,iTrial)),kMax,rmseMax); %original/full data
                    % [higBL1(iSide,iElec,iTrial,1), higBL1(iSide,iElec,iTrial,2), higBL1(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBl(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                    % hurst1(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBl(iSide,iElec,:,iTrial)));
                    % for iFilt = 1:length(bandPassFilterFreq)

                    if fastFlag
                        [higBL2{iFilt}(iSide,iElec,iTrial,1),~,~,lnK{iFilt}(iSide,iElec,iTrial,:),lnL{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);
                    else
                        [higBL2{iFilt}(iSide,iElec,iTrial,1),higBL2{iFilt}(iSide,iElec,iTrial,2),higBL2{iFilt}(iSide,iElec,iTrial,3),lnK{iFilt}(iSide,iElec,iTrial,:),lnL{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);
                    end
                    %[higBL2{iFilt}(iSide,iElec,iTrial,1),higBL2{iFilt}(iSide,iElec,iTrial,2),higBL2{iFilt}(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);

                    % hurst2{iFilt}(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)));
                    % end

                    if stFlag
                        %[higStOr(iSide,iElec,iTrial,1),higStOr(iSide,iElec,iTrial,2),higStOr(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(stData(iSide,iElec,:,iTrial)),kMax);
                        % [higSt1(iSide,iElec,iTrial,1),higSt1(iSide,iElec,iTrial,2),higSt1(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewSt(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                        % hurstSt1(iSide,iElec,iTrial) = genhurst(squeeze(dataNewSt(iSide,iElec,:,iTrial)));
                        % for iFilt = 1:length(bandPassFilterFreq)
                        if fastFlag
                            [higSt2{iFilt}(iSide,iElec,iTrial,1),~,~,lnKSt{iFilt}(iSide,iElec,iTrial,:),lnLSt{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);
                        else
                            [higSt2{iFilt}(iSide,iElec,iTrial,1),higSt2{iFilt}(iSide,iElec,iTrial,2),higSt2{iFilt}(iSide,iElec,iTrial,3),lnKSt{iFilt}(iSide,iElec,iTrial,:),lnLSt{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag,optimiseRmseFlag);
                        end
                        %  hurstSt2{iFilt}(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)));
                        %  end
                    end
                end

            else
                blData(iSide,iElec,1:numTrials,1:numTimes) = nan;
                % dataNewBl(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                %higBLOr(iSide,iElec,[1:numTrials],[1:3]) = nan;
                %higBL1(iSide,iElec,[1:numTrials],[1:3]) = nan;
                % hurst1(iSide,iElec,[1:numTrials]) = nan;
                % for iFilt = 1:length(bandPassFilterFreq)
                higBL2{iFilt}(iSide,iElec,1:numTrials,1) = nan;
                lnK{iFilt}(iSide,iElec,1:numTrials,1:kMax(iFilt)) = nan;
                lnL{iFilt}(iSide,iElec,1:numTrials,1:kMax(iFilt)) = nan;
                dataNewBandBl{iFilt}(iSide,iElec,1:numTimes,1:numTrials) = nan;
                % hurst2{iFilt}(iSide,iElec,[1:numTrials]) = nan;
                % end

                if stFlag
                    stData(iSide,iElec,1:numTrials,1:numTimes) = nan;
                    % dataNewSt(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                    %higStOr(iSide,iElec,[1:numTrials]) = nan;
                    %higSt1(iSide,iElec,[1:numTrials],[1:3]) = nan;
                    % hurstSt1(iSide,iElec,[1:numTrials]) = nan;
                    % for iFilt = 1:length(bandPassFilterFreq)
                    lnKSt{iFilt}(iSide,iElec,1:numTrials,1:kMax(iFilt)) = nan;
                    lnLSt{iFilt}(iSide,iElec,1:numTrials,1:kMax(iFilt)) = nan;
                    higSt2{iFilt}(iSide,iElec,1:numTrials,1) = nan;
                    % hurstSt2{iFilt}(iSide,iElec,[1:numTrials]) = nan;
                    dataNewBandSt{iFilt}(iSide,iElec,1:numTimes,1:numTrials) = nan;
                    % end
                end
            end
        end
        %toc;
    end
    if optimiseRmseFlag
        kMax(iFilt) = higBL2{iFilt}(iSide,iElec,iTrial,3);
        kMax(iFilt+1) = kMax(iFilt);
    end
end
% close(hW);
% for iFilt = 1:length(bandPassFilterFreq)
%     dataNewBandBl{iFilt} = permute(squeeze(dataNewBandBl{iFilt}),[1 3 2]);
%     lnK{iFilt} = squeeze(lnK{iFilt});
%     lnL{iFilt} = squeeze(lnL{iFilt});
%     higBL2{iFilt} = squeeze(higBL2{iFilt});
%     hurst2{iFilt} = squeeze(hurst2{iFilt});
%
%     if stFlag
%         dataNewBandSt{iFilt} = permute(squeeze(dataNewBandSt{iFilt}),[1 3 2]);
%         lnKSt{iFilt} = squeeze(lnKSt{iFilt});
%         lnLSt{iFilt} = squeeze(lnLSt{iFilt});
%         higSt2{iFilt} = squeeze(higSt2{iFilt});
%         hurstSt2{iFilt} = squeeze(hurstSt2{iFilt});
%     end
% end
for iFilt = 1:length(bandPassFilterFreq)
    dataNewBandBl{iFilt} = squeeze(nanmean(dataNewBandBl{iFilt},4));
    lnK{iFilt} = squeeze(nanmean(lnK{iFilt},3));
    lnL{iFilt} = squeeze(nanmean(lnL{iFilt},3));
    higBL2{iFilt} = squeeze(nanmean(higBL2{iFilt},3));
    % hurst2{iFilt} = squeeze(hurst2{iFilt});

    if stFlag
        dataNewBandSt{iFilt} = squeeze(nanmean(dataNewBandSt{iFilt},4));
        lnKSt{iFilt} = squeeze(nanmean(lnKSt{iFilt},3));
        lnLSt{iFilt} = squeeze(nanmean(lnLSt{iFilt},3));
        higSt2{iFilt} = squeeze(nanmean(higSt2{iFilt},3));
        %hurstSt2{iFilt} = squeeze(hurstSt2{iFilt});
    end
end
data.blData = squeeze(blData);
%data.dataNewBl = permute(squeeze(dataNewBl),[1 3 2]);
data.dataNewBandBl = dataNewBandBl;

%fracDim.higBLOr = squeeze(higBLOr);
%fracDim.higBL1 = squeeze(higBL1);
fracDim.kmax = kMax;
fracDim.minusLnK = lnK;
fracDim.lnL = lnL;
fracDim.higBL2 = higBL2;
fracDim.bandPassFilterFreq = bandPassFilterFreq;
fracDim.centreFreq = centreFreq;
%fracDim.hurst1 = squeeze(hurst1);
%fracDim.hurst2 = hurst2;

%fracDim.bandPassFilterFreqHurst = bandPassFilterFreq;

if stFlag
    data.stData = squeeze(stData);
    % data.dataNewSt = permute(squeeze(dataNewSt),[1 3 2]);
    data.dataNewBandSt = dataNewBandSt;
    %fracDim.higStOr = squeeze(higStOr);
    %fracDim.higSt1 = squeeze(higSt1);
    fracDim.minusLnKSt = lnKSt;
    fracDim.lnLSt = lnLSt;
    fracDim.higSt2 = higSt2;
    %fracDim.hurstSt1 = squeeze(hurstSt1);
    %fracDim.hurstSt2 = hurstSt2;

end
end