% analyseAndSaveValuesIndividualSubject
% This program calculates and saves PSDs and TF plots that are needed to
% generate all plots 
% Call this function from runAnalyseAndSaveValuesIndividualSubject.m

function val= getDataForAllSubject(folderSourceString,folderOutString, subjectName,projectName,subProjectName,refType,protocolType,parameter,cleanDataFlag)
 if ~exist('stRange','var');         stRange = [0.25 0.75];              end
if ~exist('removeMicroSaccadesFlag','var'); removeMicroSaccadesFlag=0;  end
if ~exist('cleanDataFlag','var')    cleanDataFlag = 0;  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(protocolType,'TFCP')
    condVals=16; % SSEVEPFreq/2
else
    condVals=[];
end

if cleanDataFlag
     dataFolder = fullfile(folderSourceString,'cleanData',protocolType); 
else
    dataFolder = fullfile(folderSourceString,'decimatedData',projectName,protocolType);
end

% if strcmp(protocolType,'eyes_closed') || strcmp(protocolType,'eyes_open') || strcmp(protocolType,'eyes_closed_v5')
%     stFlag=0;   %existence of stimulus
%     analyzedDataFolder = fullfile(folderOutString,['analyzedData' protocolType],projectName,protocolType);
% else
%     stFlag=1;
%     analyzedDataFolder = fullfile(folderOutString,'analyzedDataFromCleanData',projectName,protocolType);
% end

if strcmp(protocolType,'eyes_closed') || strcmp(protocolType,'eyes_open') || strcmp(protocolType,'eyes_closed_v5')
    stFlag=0;   %existence of stimulus
    analyzedDataFolder = fullfile(folderOutString,['analyzedData' protocolType],projectName,subProjectName);
else
    stFlag=1;
    analyzedDataFolder = fullfile(folderOutString,'analyzedDataFromCleanData',projectName,subProjectName);
end

%analyzedDataFolder = fullfile(pwd,['analyzedData' protocolType],projectName,protocolType); % analysedFolder now in local project directory
makeDirectory(analyzedDataFolder);

% Analysis file
% if removeMicroSaccadesFlag
%     analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
%         '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_NoMS.mat']);
%     % numMSInRangePerProtocol ; % TODO
% else
%     analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
%         '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '.mat']);
%     numMSInRangePerProtocol = [];
% end

%analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType 'withTimeSeries.mat']);
analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '.mat']);

if isfile(analysisDetailsFile)
x = load(analysisDetailsFile,parameter);
val = x.(parameter);
else
    disp('analysis file does not exist')
    val = {nan(64,1)};
    assignin('base',parameter,val)

end

%     numMSInRangePerProtocol = [];
% 
% 
% [expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);
% 
% if usableDataFlag && ~isempty(expDates)
%     clear fileLists
%     for iProt = 1:length(expDates)
%         fileLists{iProt} = [subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat']; %#ok<AGROW>
%         badTrialFileLists{iProt} = [subjectName '_' protocolNames{iProt} '_badTrials_DD.mat'];
%     end
    
    % 1. Calculate PSDs and TF plots for selected electrodes
           %  electrodeList = getElectrodeList(capLayout{1},refType); %selected electrodes
    
%     [allProtocolsBLData,stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erpData,timeVals,numGoodTrials,numAnalysedElecs]=...
%         getDataSingleSubject(dataFolder,fileLists,capLayout,electrodeList,stRange,1,numMSInRangePerProtocol,condVals);

 %      save(analysisDetailsFile,'allProtocolsBLData','stPowerVsFreq','blPowerVsFreq',...
%         'freqVals','tfPower','timeValsTF','freqValsTF','erpData','timeVals',...
%         'numGoodTrials','numAnalysedElecs','allProtocolsBLDataTopo','stPowerVsFreqTopo','blPowerVsFreqTopo','freqValsTopo');


%     
    % 2. Calculate PSDs for topoplots - TF data is not stored
    
%          electrodeList = getElectrodeList(capLayout{1},refType,1);% all electrodes

         %% for fooof analysis
%     [allProtocolsBLDataTopo,stPowerVsFreqTopo,blPowerVsFreqTopo,freqVals,~,~,~,~,timeVals,numGoodTrials,numAnalysedElecs]=...
%         getDataSingleSubject_V2(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);

% save(analysisDetailsFile,...
%         'timeVals', 'numGoodTrials','numAnalysedElecs','allProtocolsBLDataTopo','stPowerVsFreqTopo','blPowerVsFreqTopo','freqVals');
% 

%% for fracDim analysis

%[~,fracDim] = getDataSingleSubject_FracDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);
%save(analysisDetailsFile, 'fracDim');

%lag = getDataSingleSubject_CorrDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);
%save(analysisDetailsFile, 'lag');
% dim = getDataSingleSubject_CorrDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);
% save(analysisDetailsFile, 'dim','-append');

%    Not saving time series right now
%    

end