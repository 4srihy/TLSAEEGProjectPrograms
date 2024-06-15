% analyseAndSaveValuesIndividualSubject
% This program calculates and saves PSDs and TF plots that are needed to
% generate all plots 
% Call this function from runAnalyseAndSaveValuesIndividualSubject.m

function analyseAndSaveValuesIndividualSubject_v2(folderSourceString,folderOutString, subjectName,projectName,subProjectName,refType,protocolType,stRange,removeMicroSaccadesFlag,cleanDataFlag,freqRange,visElecFlag)
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
    if cleanDataFlag
        analyzedDataFolder = fullfile(folderOutString,'analyzedDataFromCleanData',projectName,subProjectName);
    else
          analyzedDataFolder = fullfile(folderOutString,'analyzedDataFromDecimatedData',projectName,subProjectName);
    end
end


%analyzedDataFolder = fullfile(pwd,['analyzedData' protocolType],projectName,protocolType); % analysedFolder now in local project directory
makeDirectory(analyzedDataFolder);

% Analysis file
%% PSD
% if removeMicroSaccadesFlag
%     analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
%         '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_NoMS.mat']);
%     % numMSInRangePerProtocol ; % TODO
% else
%     analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType ...
%         '_stRange_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '.mat']);
%     numMSInRangePerProtocol = [];
% end

%%  HFD
%analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType 'withTimeSeries.mat']);
switch subProjectName
    case 'higFracDim'
        %analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '.mat']);
        analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_freq_1_90_kmax250.mat']);
        %analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_freqWidth100_kmax5.mat']);
        %analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_TF_nonoverlapping.mat']);
       % analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_TF_nonoverlapping_kmax5.mat']);
    case 'CorrDim'   
        %% CD
        if visElecFlag
            analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_' num2str(freqRange(1)) 'to' num2str(freqRange(2)) 'HzVisElecs.mat']);
        else
            analysisDetailsFile = fullfile(analyzedDataFolder,[subjectName '_' refType '_' num2str(freqRange(1)) 'to' num2str(freqRange(2)) 'Hz.mat']);
        end
end
    numMSInRangePerProtocol = [];


[expDates,protocolNames,capLayout,usableDataFlag] = getProtocolDetailsForAnalysis(projectName,subjectName,protocolType);

if usableDataFlag && ~isempty(expDates)
    clear fileLists
    for iProt = 1:length(expDates)
        fileLists{iProt} = [subjectName '-' expDates{iProt} '-' protocolNames{iProt} '.mat']; %#ok<AGROW>
        badTrialFileLists{iProt} = [subjectName '_' protocolNames{iProt} '_badTrials_DD.mat'];
    end
end   
    % 1. Calculate PSDs and TF plots for selected electrodes
           %  electrodeList = getElectrodeList(capLayout{1},refType); %selected electrodes
    
%     [allProtocolsBLData,stPowerVsFreq,blPowerVsFreq,freqVals,tfPower,timeValsTF,freqValsTF,erpData,timeVals,numGoodTrials,numAnalysedElecs]=...
%         getDataSingleSubject(dataFolder,fileLists,capLayout,electrodeList,stRange,1,numMSInRangePerProtocol,condVals);

 %      save(analysisDetailsFile,'allProtocolsBLData','stPowerVsFreq','blPowerVsFreq',...
%         'freqVals','tfPower','timeValsTF','freqValsTF','erpData','timeVals',...
%         'numGoodTrials','numAnalysedElecs','allProtocolsBLDataTopo','stPowerVsFreqTopo','blPowerVsFreqTopo','freqValsTopo');


%     
    % 2. Calculate PSDs for topoplots - TF data is not stored
    
         electrodeList = getElectrodeList(capLayout{1},refType,1);% all electrodes

         %% for fooof analysis
%     [allProtocolsBLDataTopo,stPowerVsFreqTopo,blPowerVsFreqTopo,freqVals,~,~,~,~,timeVals,numGoodTrials,numAnalysedElecs]=...
%         getDataSingleSubject_V2(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);

% save(analysisDetailsFile,...
%         'timeVals', 'numGoodTrials','numAnalysedElecs','allProtocolsBLDataTopo','stPowerVsFreqTopo','blPowerVsFreqTopo','freqVals');
% 

%% for Higuchi fracDim analysis
tic;
[~,fracDim] = getDataSingleSubject_FracDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);
toc;
if isfile(analysisDetailsFile)
    l = load(analysisDetailsFile);
    fd = struct2cell(l.fracDim);
    fdName = fieldnames(l.fracDim);

    fd_new = struct2cell(fracDim);
    fd_newName = fieldnames(fracDim);

    for iField = 1:length(fd_new)
        if ~isempty(find(strcmp(fdName,fd_newName{iField})==1))
            if iscell(fd_new{iField})
                fieldOldPos = strcmp(fdName,fd_newName{iField});
                fd{fieldOldPos} = [fd{fieldOldPos}, fd_new{iField}];
            end
        else
            fd = [fd; fd_new{iField}];
            fdName = [fdName; fd_newName{iField}];
    end
    fracDim = cell2struct(fd,fdName);
end
end
save(analysisDetailsFile, 'fracDim');
%% for corrdim
%lag = getDataSingleSubject_CorrDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType);
%save(analysisDetailsFile, 'lag');
%fracDim = getDataSingleSubject_CorrDim(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType,freqRange);

% fracDim = getDataSingleSubject_CorrDimStBl(dataFolder,fileLists,badTrialFileLists,capLayout,electrodeList,stRange,0,numMSInRangePerProtocol,condVals,stFlag,refType,freqRange,visElecFlag);
% save(analysisDetailsFile, 'fracDim');


%    Not saving time series right now
%    


end