%by Srishty 12-07-22

function [stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erp,numAnalysedElecs]=getDataSingleSessionWithTrials(eegData,timeVals,electrodeList,stRange,TFFlag,params,stFlag,refType)

if stFlag
    blRange = [-diff(stRange) 0];  
else
    blRange = [timeVals(1)  timeVals(length(timeVals))];     stRange = blRange;
end

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

% Set MT parameters
if ~exist('params','var') || isempty(params)
    params.tapers   = [1 1];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 1000];
    params.trialave = 0;
end

movingwin = [0.25 0.025];

% Initialize
numElectrodes = length(electrodeList{1});
numSides = length(electrodeList);

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


for iElec=1:numElectrodes % For each electrode or electrode pair
    for iSide=1:numSides    % For each side
        
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
        erp(iSide,iElec,:) = mean((eeg - repmat(mean(eeg(:,blPos),2),1,size(eeg,2))),1); % Correct for DC Shift (baseline correction)
        
        if TFFlag
            params.trialave=1;
            [tfPower(iSide,iElec,:,:),timeValsTF0,freqValsTF] = mtspecgramc(eeg',movingwin,params);
            timeValsTF = timeValsTF0 + timeVals(1);
        else
            tfPower = [];
            timeValsTF = [];
            freqValsTF = [];
        end
        if stFlag   stPowerVsFreq(iSide,iElec,:,:)= mtspectrumc(squeeze(eeg(:,stPos))',params); end
        blPowerVsFreq(iSide,iElec,:,:)= mtspectrumc(squeeze(eeg(:,blPos))',params);
        
        if all(isnan(eeg(:))) % Discard bad electrodes
            erp(iSide,iElec,:) = NaN(1,size(erp,3));
            if TFFlag; tfPower(iSide,iElec,:,:) = NaN(size(tfPower,3),size(tfPower,4)); end
             if stFlag stPowerVsFreq(iSide,iElec,:,:) = NaN(1,size(stPowerVsFreq,3),size(stPowerVsFreq,4)); end
            blPowerVsFreq(iSide,iElec,:,:) = NaN(1,size(blPowerVsFreq,3),size(blPowerVsFreq,4));
            numAnalysedElecs(iSide) = numAnalysedElecs(iSide)-1;
        end
    end
end
% close(hW);

% Find mean across electrodes
erp = removeDimIfSingleton(nanmean(erp,2),2);
tfPower = removeDimIfSingleton(nanmean(tfPower,2),2);
stPowerVsFreq = removeDimIfSingleton(nanmean(stPowerVsFreq,2),2);
blPowerVsFreq = removeDimIfSingleton(nanmean(blPowerVsFreq,2),2);

if ~isempty(eeg)
    [~,freqVals]= mtspectrumc(squeeze(eeg(1,stPos))',params);
else
    [~,freqVals]= mtspectrumc(squeeze(eeg(:,stPos))',params);
end

end