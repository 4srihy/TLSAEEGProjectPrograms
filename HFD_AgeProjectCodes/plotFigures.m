clear;
folderSourceString{1} ='C:\Users\Srishty\OneDrive - Indian Institute of Science\Data\TLSA\analyzedData'; % Indicate the folder of analyzedData
figNum = 1; %figure no. to be plotted

projectName = 'ADGammaProject'; % Only this dataset, which is the main TLSA dataset, is configured as of now. Other options - 'AgeProjectRound1' and 'VisualGamma' may not work
stRange = [0.25 0.75];

% Choose one of these options
refType = 'unipolar'; % 'unipolar' % Set reference type here.
protocolType = 'SF_ORI'; % 'TFCP'; % SF_ORI for gamma, TFCP for SSVEP, 'eyes_closed'
removeMicroSaccadesFlag = 0; % 0 or 1.
saveFileFlag = 0;

goodSubjectsList = getGoodSubjectsProjectwise(projectName);
uniqueSubjectNames = getGoodFileNamesForSubjects(goodSubjectsList{1});
[~,~,cdrList] = getDemographicDetails(projectName,uniqueSubjectNames);

%ageLim = 5;
healthyPos = strcmp(cdrList,'HV');

goodSubjectsAll = uniqueSubjectNames(healthyPos);

goodSubjectsAll(99) = []; %remooving 275KK
%removing subjects with captype ActicapPosterioir31
detail = load("ADGammaProjectDetails.mat");
rejectCapPos=find(strcmp(detail.expDetails(:,3),'actiCap31Posterior'));
subNames = detail.expDetails(rejectCapPos,1);
rejectSubPos= [];
for i=1:length(subNames)
    k = find(strcmp(goodSubjectsAll,subNames(i)));
    rejectSubPos = [rejectSubPos k];
end
goodSubjectsAll(rejectSubPos) = [];

[ageList,genderList] = getDemographicDetails(projectName,goodSubjectsAll);

ageLims = [50 65 90];
clear subjectNameList
for iAge = 1:(length(ageLims)-1)
    ageGroupPos{iAge} = (ageList>=ageLims(iAge) & (ageList<ageLims(iAge+1)));
    subjectNameList{iAge} = goodSubjectsAll(ageGroupPos{iAge});
    %strList{iAge} = ['(' num2str(ageLims(iAge)) '-' num2str(ageLims(iAge+1)) ')'];
end
strList{1} = 'Mid'; strList{2} = 'Old';

%%%%%%%%%%%%%%%%%%%% OPTIONAL: CAN DO GENDERWISE ANALYSIS %%%%%%%%%%%%%%%
%by giving subjectNameList as subjectNameListGender
%or by giving subjectNameList as subMale as\nd subFemale for young/old group

% subjectNameListGender{1}=[];    subjectNameListGender{2}=[];
%  for i = 1:(length(ageLims)-1)
%             [ageListnew{i}, genderListnew{i}] = getDemographicDetails(projectName, subjectNameList{i});
%              MalePos{i} = strcmp('M',genderListnew{i});
%             subMale{i} = subjectNameList{i}(MalePos{i});
%              FemalePos{i} = strcmp('F',genderListnew{i});
%             subFemale{i} = subjectNameList{i}(FemalePos{i});
%
%             subjectNameListGender{1} = [subjectNameListGender{1} subMale{i}];
%             subjectNameListGender{2} = [subjectNameListGender{2} subFemale{i}];
%    end

%subjectNameList = subjectNameListGender;%subFemale;%
%subjectNameList{1} = subMale{2};    subjectNameList{2} = subFemale{2};
%strList{1} = 'M Old'; strList{2} = 'F Old';
%strList{1} = 'M'; strList{2} = 'F';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% running the code for figure
if figNum<8
    run(['fig' num2str(figNum) '.m']);
else
    run('suppFig1.m');
end
