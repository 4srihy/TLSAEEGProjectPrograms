%by srishty 30-05-23
function [data,fracDim,timeVals,numGoodTrials,numAnalysedElecs]=getDataSingleSubject_FracDim(cleanDataFolder,fileLists,badTrialFileLists,capType,electrodeList,stRange,TFFlag,numMSInRangePerProtocol,condVals,stFlag,refType,discardBadElecFlag)

if ~exist('TFFlag','var') || isempty(TFFlag); TFFlag= 1; end
if ~exist('stFlag','var');       stFlag = 1;     end
if ~exist('stRange','var') || isempty(stRange); stRange = [0.25 0.75]; end
if ~exist('params','var'); params=[]; end

if exist('numMSInRangePerProtocol','var') && ~isempty(numMSInRangePerProtocol)
    if length(numMSInRangePerProtocol)~=length(fileLists); error('Remove bad protocols first...'); end
end

if stFlag
    rmsBadInfoFolder = fullfile(pwd,'badTrialsForDecimatedDataIncludingRMS');
else
    rmsBadInfoFolder = fullfile(pwd,'badTrialsForEyesClosedIncludingRMS'); %05-11-22
end

iLoop = 1;%
for iProt = 1:length(fileLists)
    sessionData = load(fullfile(cleanDataFolder,fileLists{iProt}),'eegData','badElecs','trialConditionVals','timeVals');
    eegData = sessionData.eegData;
    badElecs = sessionData.badElecs;
    trialConditionVals = sessionData.trialConditionVals;
    timeVals = sessionData.timeVals;
    if ~isempty(trialConditionVals) && length(trialConditionVals)~=size(eegData,2)
        error('Wrong number of trials in EEG data or condition values...');
    end

    % removing rms bad trials..only for clean data and decimated data for SF_ORI
     %05-11-22
        rmsBadElecsNTrials = load(fullfile(rmsBadInfoFolder,badTrialFileLists{iProt}),'badElecs','badTrials');
        badElecs = rmsBadElecsNTrials.badElecs;
        badTrials = rmsBadElecsNTrials.badTrials;
        eegData(:,badTrials,:) =[];
    
    % Check for bad protocols
    badImpedanceElecs = badElecs.badImpedanceElecs;
    noisyElecs = badElecs.noisyElecs;
    flatPSDElecs = badElecs.flatPSDElecs;
    allBadElecs = (union(union(badImpedanceElecs,noisyElecs),flatPSDElecs))';

    electrodeListVis = getElectrodeList(capType,refType);
    clear goodSideFlag
    for iSide = 1:length(electrodeListVis)
        clear goodElecFlag
        for iBipElec = 1:length(electrodeListVis{iSide})
            goodElecFlag(iBipElec) = ~any(ismember(electrodeListVis{iSide}{iBipElec},allBadElecs)); %#ok<*AGROW>
        end
        if any(goodElecFlag)
            goodSideFlag(iSide) = true;
        else
            goodSideFlag(iSide) = false;
        end
    end

    % Deal with bad electrodes: discard data unless specified
    if exist('discardBadElecFlag','var') && ~discardBadElecFlag
        warning('Not discarding data from bad electrodes...')
    else
        for iBE = 1:length(allBadElecs)
            eegData(allBadElecs(iBE),:,:) = nan;%NaN(size(eegData,2),size(eegData,3));
        end
    end

    % Remove microsaccde containing trials from analysis
    if exist('numMSInRangePerProtocol','var') && ~isempty(numMSInRangePerProtocol)
        msTrials = numMSInRangePerProtocol{iProt}>0;
        eegData(:,msTrials,:) = [];
        trialConditionVals(msTrials) = [];
    end

    if exist('condVals','var') && ~isempty(condVals); eegData(:,~ismember(trialConditionVals,condVals),:) = []; end
    if all(goodSideFlag)
        goodProtFlag(iLoop)=true;
    else
        goodProtFlag(iLoop)=false;
    end

    %%%%%%%%%%%%
    %[stPowerVsFreq(iLoop,:,:),blPowerVsFreq(iLoop,:,:),freqVals,tfPower(iLoop,:,:,:),timeValsTF,freqValsTF,erp(iLoop,:,:,:),numAnalysedElecs(iLoop,:)]= getDataSingleSession(eegData,timeVals,electrodeList,stRange,TFFlag,params);
    %     [stPowerVsFreq0{iLoop},blPowerVsFreq0{iLoop},freqVals,tfPower(iLoop,:,:,:),timeValsTF,freqValsTF,erp(iLoop,:,:,:),numAnalysedElecs(iLoop,:)]=...
    %         getDataSingleSessionWithTrials(eegData,timeVals,electrodeList,stRange,TFFlag,params,stFlag,refType); %by Srishty 12-07-22
    [~,fracDim0{iLoop},numAnalysedElecs] = getHigFracDimSingleSessionWithTrials(eegData,timeVals,electrodeList,stRange,stFlag,refType);

    %[data0{iLoop},fracDim0{iLoop},numAnalysedElecs] = getHigFracDimSingleSessionWithTrialsContinuousFreqRange(eegData,timeVals,electrodeList,stRange,stFlag,refType);

    %%%%%%%% for tf data %%%%%%%%%%%%%
    %         electrodeList = 1:64;    data = [];
    %     [fracDim0{iLoop}] = genTFPlotForHFD(eegData,timeVals,electrodeList,refType,100,12,250);

    %%%%%%%%%%%%%%%%%%
    numGoodTrials(iLoop) = size(eegData,2);
    iLoop = iLoop+1;
end

%

if any(goodProtFlag)
    catDim = 3;  %check which one is better 1 or 2; 1 for genTFandOtherFreqWidth with mean trial data, 3 for rest
    data = [];
    if iLoop>2  %need to change ..can't run for baseline and stimulus
        fd1 = struct2cell(fracDim0{1});
        fdName = fieldnames(fracDim0{1});

        if exist('data0','var')
            d1 = struct2cell(data0{1});
            dName = fieldnames(data0{1});
        end
        for i=1:iLoop-2
            fd2 = struct2cell(fracDim0{i+1});
            for j=1:length(fd1)
                if iscell(fd1{j})
                    for k = 1:length(fd1{j})
                        fd1{j}{k} = cat(catDim,fd1{j}{k},fd2{j}{k});
                    end
                else
                    fd1{j} = cat(catDim,fd1{j},fd2{j});
                end
            end

            if exist('data0','var')
                d2 = struct2cell(data0{i+1});
                for j=1:length(d1)
                    if iscell(d1{j})
                        for k = 1:length(d1{j})
                            d1{j}{k} = cat(2,d1{j}{k},d2{j}{k});
                        end
                    else
                        d1{j} = cat(2,d1{j},d2{j});
                    end
                end
            end
        end
        if exist('data0','var')   data = cell2struct(d1,dName);   end
        fracDim = cell2struct(fd1,fdName);

    else
        if exist('data0','var')    data = data0{1}; end
        fracDim = fracDim0{1};
    end
    %     stPowerVsFreq = removeDimIfSingleton(mean(stPowerVsFreq(goodProtFlag,:,:),1),1);
    %     blPowerVsFreq = removeDimIfSingleton(mean(blPowerVsFreq(goodProtFlag,:,:),1),1);

else
    disp('bad visual electrodes: data not saved');
    data = [];
    fracDim = [];

end
end

%by Srishty 12-07-22

