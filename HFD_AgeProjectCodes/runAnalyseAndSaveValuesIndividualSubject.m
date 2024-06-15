% Modified from the program written by Murty Dinavahi (MD)

%further edited by Srishty 23-8-22
%now can work on decimated data as well 
% (need to change the details in 'analyzeAndSaveIndividualSubjects')
clear;

% Mandatory fixed options

folderSourceString = '\\NEOLABDATA\NeoLabData\Projects\TLSAEEGProject'; % Source of clean data
folderOutString =  'G:\OneDrive - Indian Institute of Science\Data\TLSA'; %Path for saved data
projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];
%freqRange = [7 13];

% Choose one of these options
refType = 'unipolar'; % 'unipolar', 'bipolar or 'average' % Set reference type here.               
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP,eyes_closed
subFolderName =  "higFracDim";%'CorrDim'; %
freqRange = [80,140];
visElecFlag = 0;  %1 : only visual electrodes, 0: all electrodes

cleanDataFlag = 1; %0 for decimated data ,1 for clean data
removeMicroSaccadesFlag = 0; % 0 or 1.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subProjectName = 'age';
goodSubjects = getGoodSubjectsProjectwise(subProjectName,1);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjects{1});

subjectswithIssueEyesclosed = [152 154 159 161 165 166 168 169 176];
subjectswithIssueEyesOpen = [19 subjectswithIssueEyesclosed ];
parameter = 'dim';
for iSub = setdiff(121:150,subjectswithIssueEyesclosed)%176:200%length(uniqueSubjectNames)%setdiff(1:length(uniqueSubjectNames),subjectswithIssueEyesclosed)
    subjectName = uniqueSubjectNames{iSub};
    disp([num2str(iSub) ': ' subjectName]);
    analyseAndSaveValuesIndividualSubject_v2(folderSourceString, folderOutString, subjectName,projectName,subFolderName,refType,protocolType,stRange,removeMicroSaccadesFlag,cleanDataFlag,freqRange,visElecFlag); % Save data in analyzedData
end

 