function [data,fracDim,numAnalysedElecs]=getHigFracDimSingleSessionWithTrials_includingHurst(eegData,timeVals,electrodeList,stRange,stFlag,refType)

%if ~exist('kMax','var')         kMax  = 20;         end
bandPassFilterFreq = {[1 90]};%{[35 66],[210 290]};% {[8 12], [20 35],[60 90], [110 140], [160 190],[210 240],[260 290]};
kMax = [20];%[5 5];%kmax should have size of length(bandpassFilterFreq)%[10];%[10 10 10 7 5 4 3]; % for 1-90 :kmax=10
%bandPassFilterFreq = {[7 13], [20 35],[60 90]}; % alpha, slow gamma, and aperiodic
awayFromCentre = 2;
for freqAvoidcentre = 50:50:250
    stopFreqRange{freqAvoidcentre/50} = [freqAvoidcentre-awayFromCentre freqAvoidcentre+awayFromCentre];
end 

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
d{1} = designfilt('highpassiir','FilterOrder',8, 'PassbandFrequency',1,'PassbandRipple',0.2,'SampleRate',2500);
for iFilt = 1:length(bandPassFilterFreq)
    b{iFilt} = designfilt('bandpassiir','FilterOrder',4,'PassbandFrequency1',bandPassFilterFreq{iFilt}(1),'PassbandFrequency2',bandPassFilterFreq{iFilt}(2),'PassbandRipple',0.2,'SampleRate',2500);
    iStop0 = 0;
    for iStop = 1:length(stopFreqRange)
        if (stopFreqRange{iStop}(1)<bandPassFilterFreq{iFilt}(2)  && stopFreqRange{iStop}(1)>bandPassFilterFreq{iFilt}(1)) || (stopFreqRange{iStop}(2)<bandPassFilterFreq{iFilt}(2)  && stopFreqRange{iStop}(2)>bandPassFilterFreq{iFilt}(1))
            iStop0 = iStop0+1;
            bStop{iFilt}{iStop0} = designfilt('bandstopiir','FilterOrder',4,'PassbandFrequency1',stopFreqRange{iStop}(1),'PassbandFrequency2',stopFreqRange{iStop}(2),'PassbandRipple',0.2,'SampleRate',2500);
        end
    end
end
tic;
for iElec=1:numElectrodes % For each electrode or electrode pair
    
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
            dataNewBl(iSide,iElec,:,:) = filtfilt(d{1},squeeze(blData(iSide,iElec,:,:))');
            for iFilt = 1:length(bandPassFilterFreq)
                dataNewBandBl{iFilt}(iSide,iElec,:,:)   = filtfilt(b{iFilt},squeeze(blData(iSide,iElec,:,:))');
                if exist('bStop','var')
                    if ~isempty(bStop{iFilt})
                    for iStopFilt = 1:length(bStop{iFilt})
                        dataNewBandBl{iFilt}(iSide,iElec,:,:) = filtfilt(bStop{iFilt}{iStopFilt},squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,:)));
                    end
                    end
                end
            end

            if stFlag
                stData(iSide,iElec,:,:) =  eeg(:,stPos);
                dataNewSt(iSide,iElec,:,:) = filtfilt(d{1},squeeze(stData(iSide,iElec,:,:))');
                for iFilt = 1:length(bandPassFilterFreq)
                    dataNewBandSt{iFilt}(iSide,iElec,:,:)   = filtfilt(b{iFilt},squeeze(stData(iSide,iElec,:,:))');
                    if exist('bStop','var')
                    if ~isempty(bStop{iFilt})
                    for iStopFilt = 1:length(bStop{iFilt})
                        dataNewBandSt{iFilt}(iSide,iElec,:,:) = filtfilt(bStop{iFilt}{iStopFilt},squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,:)));
                    end
                    end
                end
                end
            end

            for iTrial = 1:numTrials
              %  [higBLOr(iSide,iElec,iTrial,1),higBLOr(iSide,iElec,iTrial,2),higBLOr(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(blData(iSide,iElec,:,iTrial)),kMax,rmseMax); %original/full data
                if fastFlag
                    [higBL1(iSide,iElec,iTrial,1), ~, ~,~,lnL1(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBl(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                else
                    [higBL1(iSide,iElec,iTrial,1), higBL1(iSide,iElec,iTrial,2), higBL1(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBl(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                end
                hurst1(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBl(iSide,iElec,:,iTrial)));
                for iFilt = 1:length(bandPassFilterFreq)
                    if fastFlag
                        [higBL2{iFilt}(iSide,iElec,iTrial,1),~,~,~,lnL{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag);
                    else
                        [higBL2{iFilt}(iSide,iElec,iTrial,1),higBL2{iFilt}(iSide,iElec,iTrial,2),higBL2{iFilt}(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag);
                    end
                    hurst2{iFilt}(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBandBl{iFilt}(iSide,iElec,:,iTrial)));
                end

                if stFlag
                    [higStOr(iSide,iElec,iTrial,1),higStOr(iSide,iElec,iTrial,2),higStOr(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(stData(iSide,iElec,:,iTrial)),kMax);
                    if fastFlag
                        [higSt1(iSide,iElec,iTrial,1),~,~,~,lnLSt1(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewSt(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                    else
                        [higSt1(iSide,iElec,iTrial,1),higSt1(iSide,iElec,iTrial,2),higSt1(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewSt(iSide,iElec,:,iTrial)),kMax(1),rmseMax,fastFlag);
                    end
                    hurstSt1(iSide,iElec,iTrial) = genhurst(squeeze(dataNewSt(iSide,iElec,:,iTrial)));
                    for iFilt = 1:length(bandPassFilterFreq)
                        if fastFlag
                            [higSt2{iFilt}(iSide,iElec,iTrial,1),~,~,~,lnLSt{iFilt}(iSide,iElec,iTrial,:)] = HigFracDimV2(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag);
                        else
                            [higSt2{iFilt}(iSide,iElec,iTrial,1),higSt2{iFilt}(iSide,iElec,iTrial,2),higSt2{iFilt}(iSide,iElec,iTrial,3)] = HigFracDimV2(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)),kMax(iFilt),rmseMax,fastFlag);
                        end
                        hurstSt2{iFilt}(iSide,iElec,iTrial) = genhurst(squeeze(dataNewBandSt{iFilt}(iSide,iElec,:,iTrial)));
                    end
                end
            end

        else
%             blData(iSide,iElec,[1:numTrials],[1:numTimes]) = nan;
%             dataNewBl(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
 %           higBLOr(iSide,iElec,[1:numTrials],[1:3]) = nan;
            if fastFlag
                higBL1(iSide,iElec,[1:numTrials],1) = nan;
                lnL1(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;
            else
                higBL1(iSide,iElec,[1:numTrials],[1:3]) = nan;
            end
            hurst1(iSide,iElec,[1:numTrials]) = nan;
     %       lnK1(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;

            for iFilt = 1:length(bandPassFilterFreq)
                if fastFlag
                    higBL2{iFilt}(iSide,iElec,[1:numTrials],1) = nan;
                    lnL{iFilt}(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;
                else
                    higBL2{iFilt}(iSide,iElec,[1:numTrials],[1:3]) = nan;
                end
                dataNewBandBl{iFilt}(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                hurst2{iFilt}(iSide,iElec,[1:numTrials]) = nan;
                %lnK{iFilt}(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;

            end

            if stFlag
%                 stData(iSide,iElec,[1:numTrials],[1:numTimes]) = nan;
%                 dataNewSt(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
  %              higStOr(iSide,iElec,[1:numTrials]) = nan;
                if fastFlag
                    higSt1(iSide,iElec,[1:numTrials],1) = nan;
                    lnLSt1(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;
                else
                    higSt1(iSide,iElec,[1:numTrials],[1:3]) = nan;
                end
                hurstSt1(iSide,iElec,[1:numTrials]) = nan;
                lnKSt1(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;

                for iFilt = 1:length(bandPassFilterFreq)
                    if fastFlag
                        higSt2{iFilt}(iSide,iElec,[1:numTrials],1) = nan;
                        lnLSt{iFilt}(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;
                    else
                        higSt2{iFilt}(iSide,iElec,[1:numTrials],[1:3]) = nan;
                    end
                    hurstSt2{iFilt}(iSide,iElec,[1:numTrials]) = nan;
                    dataNewBandSt{iFilt}(iSide,iElec,[1:numTimes],[1:numTrials]) = nan;
                    %lnKSt{iFilt}(iSide,iElec,[1:numTrials],[1:kMax(iFilt)]) = nan;

                end
            end
        end
    end
    
end
toc; 
% close(hW);
% for iFilt = 1:length(bandPassFilterFreq)
%     dataNewBandBl{iFilt} = permute(squeeze(dataNewBandBl{iFilt}),[1 3 2]);
%     higBL2{iFilt} = squeeze(higBL2{iFilt});
%     %hurst2{iFilt} = squeeze(hurst2{iFilt});
%     %lnK{iFilt} = squeeze(lnK{iFilt});
%     lnL{iFilt} = squeeze(lnL{iFilt});
%     if stFlag
%         dataNewBandSt{iFilt} = permute(squeeze(dataNewBandSt{iFilt}),[1 3 2]);
%         higSt2{iFilt} = squeeze(higSt2{iFilt});
%         %hurstSt2{iFilt} = squeeze(hurstSt2{iFilt});
%        % lnKSt{iFilt} = squeeze(lnKSt{iFilt});
%         lnLSt{iFilt} = squeeze(lnLSt{iFilt});
%     end
% end

%% averaging across trials

%dataNewBl = squeeze(nanmean(dataNewBl,4));

%lnK{iFilt} = squeeze(nanmean(lnK{iFilt},3));

 if fastFlag;    lnL1= squeeze(nanmean(lnL1-lnL1(:,:,:,1),3));     end
higBL1 = squeeze(nanmean(higBL1,3));
hurst1 = squeeze(nanmean(hurst1,3));

if stFlag
    dataSt = squeeze(nanmean(dataNewSt,4));
    %lnKSt{iFilt} = squeeze(nanmean(lnKSt{iFilt},3));
    if fastFlag  lnLSt1 = squeeze(nanmean(lnLSt1-lnLSt1(:,:,:,1),3));    end
    higSt1 = squeeze(nanmean(higSt1,3));
    hurstSt1 = squeeze(nanmean(hurstSt1,3));
end

for iFilt = 1:length(bandPassFilterFreq)
    dataNewBandBl{iFilt} = squeeze(nanmean(dataNewBandBl{iFilt},4));
    %lnK{iFilt} = squeeze(nanmean(lnK{iFilt},3));
    if fastFlag;     lnL{iFilt} = squeeze(nanmean(lnL{iFilt}-lnL{iFilt}(:,:,:,1),3)); end
    higBL2{iFilt} = squeeze(nanmean(higBL2{iFilt},3));
    hurst2{iFilt} = squeeze(nanmean(hurst2{iFilt},3));

    if stFlag
        dataNewBandSt{iFilt} = squeeze(nanmean(dataNewBandSt{iFilt},4));
        % lnKSt{iFilt} = squeeze(nanmean(lnKSt{iFilt},3));
        if fastFlag;  lnLSt{iFilt} = squeeze(nanmean(lnLSt{iFilt}-lnLSt{iFilt}(:,:,:,1),3));  end
        higSt2{iFilt} = squeeze(nanmean(higSt2{iFilt},3));
        hurstSt2{iFilt} = squeeze(nanmean(hurstSt2{iFilt},3));
    end
end


% data.blData = squeeze(blData);
% data.dataNewBl = dataNewBl;%permute(squeeze(dataNewBl),[1 3 2]);
data.dataNewBandBl = dataNewBandBl;

%fracDim.higBLOr = squeeze(higBLOr);
fracDim.higBL1 = higBL1;
fracDim.higBL2 = higBL2;
fracDim.bandPassFilterFreq = bandPassFilterFreq;

if fastFlag
    %fracDim.minusLnK1 = squeeze(lnK1);
    fracDim.lnL1 = lnL1;
    %fracDim.minusLnK = lnK;
    fracDim.lnL = lnL;
end
% 
% fracDim.hurst1 = hurst1;
% fracDim.hurst2 = hurst2;

if stFlag
%     data.stData = squeeze(stData);
%     data.dataNewSt = dataNewSt;%permute(squeeze(dataNewSt),[1 3 2]);
    data.dataNewBandSt = dataNewBandSt;
    %fracDim.higStOr = squeeze(higStOr);
    %fracDim.minusLnKSt1 = squeeze(lnKSt1);

    if fastFlag
        fracDim.lnLSt1 = lnLSt1;
        %fracDim.minusLnKSt = lnKSt;
        fracDim.lnLSt = lnLSt;
    end

   fracDim.higSt1 = higSt1;
    fracDim.higSt2 = higSt2;
    fracDim.hurstSt1 = hurstSt1;
    fracDim.hurstSt2 = hurstSt2;

end
end