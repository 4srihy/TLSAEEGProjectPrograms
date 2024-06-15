%% to obtain power matched subject list fromnormal data files
% folderSourceString = 
% projectName = ''
function matchedSubjectNameList = getPowerMatchedSubjects(subjectNameList,folderName,matchingType,numIter,powerMatchingBinWidth)

if ~exist('matchingType','var');             matchingType = 'DiffdB';        end%'BL','ST','DiffdB','DiffOr'
if ~exist('powerMatchingBinWidth','var');    powerMatchingBinWidth = 3;      end%dB
if ~exist('numIter','var');                  numIter = 10;                   end   %no. of iterations for matched combinations

refType = 'unipolar';
freqRangeForMatching = {[8 12],[20 36],[40 66]};

numGroup= length(subjectNameList);
for i=1:numGroup
   for j=1:length(subjectNameList{i})
       subjectName = subjectNameList{i}{j};
        fileName = fullfile(folderName,[subjectName '_' refType '_stRange_250_750.mat']);
         if ~isfile(fileName)
        disp(['fileName for subject ' subjectName ' does not exist']);
    else
        x = load(fileName);
        powerBl{i}(j,:) = x.blPowerVsFreq;
        powerSt{i}(j,:) = x.stPowerVsFreq;
         end
   end
end

freqVals = x.freqVals;

switch matchingType
    case 'BL'
        powerToMatch = powerBl;
    case 'ST'
        powerToMatch = powerSt;
    case 'DiffdB'
        for i=1:numGroup
         powerToMatch{i} = 10*(log10(powerSt{i})-log10(powerBl{i}));
        end
    case 'DiffOr'
         for i=1:numGroup
        powerToMatch{i} = powerSt{i}-powerBl{i};
         end

end


  for iFreq=1:length(freqRangeForMatching) 
      freqPos = freqVals>=freqRangeForMatching{iFreq}(1) & freqVals<=freqRangeForMatching{iFreq}(2);
      for i=1:numGroup
      powerToMatchFreq{iFreq}{i} = squeeze(sum(powerToMatch{i}(:,freqPos),2));
      end
    for iTer = 1:numIter
      matchedSubjectNameList{iFreq}{iTer} = getParameterMatchedSubjectLists(subjectNameList, powerToMatchFreq{iFreq},powerMatchingBinWidth);
    end
        disp(['no. of matched subjects in freq range' num2str(freqRangeForMatching{iFreq}(1)) ' to ' num2str(freqRangeForMatching{iFreq}(2))  ' is ' num2str(length(matchedSubjectNameList{iFreq}{1}{1}))]);
  end
end