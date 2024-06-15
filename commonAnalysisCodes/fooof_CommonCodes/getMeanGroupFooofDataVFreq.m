% get results group wise
 % the resultss doesn't contain median or mean  
 %contain multiple freq_ranges
 function [groupData,totalSubjectsAnalyzed] = getMeanGroupFooofDataVFreq(subjectList,folderSourceString, projectName, refType, freqRangeWidth, powerType, OptimR_SQ, OptimExponent,medianFlag,powerFlag,str_condition,alphaFlag)

%folderSourceString ='C:\Users\Srishty\OneDrive - Indian Institute of Science\Documents\supratim\TLS\TLSAEEGProjectPrograms-master\ADGammaProjectCodes\analyzedData'; % Indicate the folder of analyzedData
%projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
%refType = 'bipolar';

%if ~exist('freq_range','var')  freq_range  = [4,44];    end
if ~exist("freqRangeWidth",'var')   freqRangeWidth = 16;        end
if ~exist('tapers')                 tapers = 1;                 end
if ~exist('powerType','var')        powerType = 'BL';           end
if ~exist('OptimR_SQ','var')        OptimR_SQ= 0.5;             end
if ~exist('OptimExponent','var')    OptimExponent= 0.1;         end
if ~exist('medianFlag','var')       medianFlag= 0;              end
if ~exist('powerFlag','var')        powerFlag= 0;               end
if ~exist('alphaFlag','var')        alphaFlag = 0;              end

%stdThresh = 3.5;
declaredBadElecs = [];
Exponent = []; Offset = []; Knee = []; specPower = []; R_SQ = [];
numSubjects = length(subjectList);
analyzedDataFolder = fullfile(folderSourceString,projectName,'FOOOF');
totalSubjectsAnalyzed = 0;
for iSub = 1:numSubjects
    subjectName = subjectList{iSub};
    %analysedFile = fullfile(analyzedDataFolder,[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '.mat']);
    %needs to be changed
   
    %analysedFile = fullfile(analyzedDataFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType 'withRandomTrials(' str_condition '.mat']);
    analysedFile = fullfile(analyzedDataFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType str_condition '.mat']);
      
   %analysedFile = fullfile(analyzedDataFolder,[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '_' powerType '.mat']);

    if ~isfile(analysedFile)
        disp(['fileName for subject ' subjectName ' does not exist']);
       
    else
        x = load(analysedFile);
        for iFreq = 1:length(x.freq_range)
            rejectElec1 = find(x.r_SQ(:,iFreq)<OptimR_SQ);
            rejectElec2 = find(x.exponent(:,iFreq)<OptimExponent);
    
            rejectElecFinal = union(union(rejectElec1,rejectElec2),declaredBadElecs);

%         logPower = log10(x.SpecPower);
%         powerMinus = nanmean(logPower,1)-stdThresh.*nanstd(logPower,1);
%         powerPlus = nanmean(logPower,1)+stdThresh.*nanstd(logPower,1);
%         numElec = size(logPower,1);
%         errElectrodes = sum((logPower(:,[5:20 30:45])<ones(numElec,1)*powerMinus(:,[5:20 30:45]) | logPower(:,[5:20 30:45])>ones(numElec,1)*powerPlus(:,[5:20 30:45])),2);
%         errElectrodeNum = find(errElectrodes~=0);
%         rejectElecFinal = union(rejectElecFinal,errElectrodeNum);

            x.exponent(rejectElecFinal,iFreq) = nan;
            x.offset(rejectElecFinal,iFreq) = nan;
            x.knee(rejectElecFinal,iFreq) = nan;
            x.r_SQ(rejectElecFinal,iFreq) = nan;
        end
      % if powerFlag     x.SpecPower(rejectElecFinal,:) = nan; end

        Exponent(iSub,:,:) = x.exponent;
        Offset(iSub,:,:) = x.offset;
        Knee(iSub,:,:) = x.knee;
        R_SQ(iSub,:,:) = x.r_SQ;
       if powerFlag     specPower(iSub,:,:) = x.SpecPower;  end
       if alphaFlag 
           [CF(iSub),PH(iSub),BW(iSub),AP(iSub),APOld(iSub)] = getAlphaParameters(x.fooof_results,x.freqVals,x.SpecPower,medianFlag);
       end
       totalSubjectsAnalyzed = totalSubjectsAnalyzed+1;
    end
end
    Exponent(find(Exponent==0)) = nan;   Offset(find(Offset==0)) = nan;     Knee(find(Knee==0)) = nan; R_SQ(find(R_SQ==0)) = nan;
    groupData.freq_range = x.freq_range;
    groupData.centreFreq  = x.centreFreq;
    groupData.freq = x.freqVals;
    groupData.exponent = Exponent;
        groupData.offset = Offset;
        groupData.knee = Knee;
        groupData.r_SQ = R_SQ;
       if powerFlag  specPower(find(specPower==0)) = nan;   groupData.SpecPower = specPower;  end
       if alphaFlag
           CF(find(CF==0)) = nan;   PH(find(PH==0)) = nan;  
           BW(find(BW==0)) = nan;   AP(find(AP==0)) = nan;
           groupData.CF =   CF;
           groupData.PH = PH;
           groupData.BW = BW;
           groupData.AP = AP;
           groupData.APOld = APOld;
       end
%     if medianFlag
%         groupData.exponent = nanmedian(Exponent,1);
%         groupData.offset = nanmedian(Offset,1);
%         groupData.knee = nanmedian(Knee,1);
%        if powerFlag  groupData.SpecPower = squeeze(nanmedian(specPower,1));  end
%     else
%         groupData.exponent = nanmean(Exponent,1);
%         groupData.offset = nanmean(Offset,1);
%         groupData.knee = nanmean(Knee,1);
%         if powerFlag  groupData.SpecPower = squeeze(nanmean(specPower,1));   end
%     end
end