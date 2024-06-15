function [freqVals,centreFreq,PSD, paramData, powerElectrodeGroup,paramElectrodeGroupData,paramSubjectElectrodeGroupData,logPSD,logPowerElectrodeGroup,totalSubjectsAnalyzed] = getGroupWiseElectrodeGroupData(subjectNameList,parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidth,tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition,alphaFlag)
    if ~exist('alphaflag','var')    alphaFlag=0; end
%folderSourceStringNew = fullfile(folderSourceString,['analyzedData' num2str(tapers)]);
    %folderSourceStringNew = fullfile(folderSourceString,'analyzedDataFromCleanData');
    %numData = getMeanGroupFooofDataVExtra(subjectNameList,folderSourceString, projectName, refType, freq_range, powerType, OptimR_SQ, OptimExponent,medianFlag,1);
    [numData,totalSubjectsAnalyzed] = getMeanGroupFooofDataVFreq(subjectNameList,folderSourceString, projectName, refType, freqRangeWidth, powerType, OptimR_SQ, OptimExponent,medianFlag,1,str_condition,alphaFlag);
    freqVals = numData.freq;
    centreFreq = numData.centreFreq;
    %Freqs = numData.freq;
    PSD = numData.SpecPower;
    paramData = numData.(parameter);

    PSD(find(PSD==0)) = nan;
    paramData(find(paramData==0)) = nan;

    for j = 1:length(electrodeGroupList)
        if medianFlag
             powerElectrodeGroup{j} = squeeze(nanmedian(PSD(:,electrodeGroupList{j},:),2));%first across electrodes for a subject and then across subjects
             paramSubjectElectrodeGroupData{j} = squeeze(nanmedian(paramData(:,electrodeGroupList{j},:),2));
            paramElectrodeGroupData{j} = nanmedian(paramSubjectElectrodeGroupData{j});
           % stdElectrodeGroupData{j} = nanstd(bootstrp(10000,@nanmedian,paramSubjectElectrodeGroupData{j}));
        else
            powerElectrodeGroup{j} =squeeze(nanmean(PSD(:,electrodeGroupList{j},:),2));
            paramSubjectElectrodeGroupData{j} = squeeze(nanmean(paramData(:,electrodeGroupList{j},:),2));
            paramElectrodeGroupData{j} = nanmean(paramSubjectElectrodeGroupData{j});
            %stdElectrodeGroupData{j} = nanstd(paramSubjectElectrodeGroupData{j})/sqrt(length(paramSubjectElectrodeGroupData{j}));
        
        end
        logPowerElectrodeGroup{j} = log10(powerElectrodeGroup{j});
    end
    
    logPSD = log10(PSD);
    
end
