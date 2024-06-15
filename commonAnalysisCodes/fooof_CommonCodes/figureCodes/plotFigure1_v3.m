%figure 1 new_v3 %1-1-23
%add additional noise
function plotFigure1_v3(subjectNameList,strList,folderSourceString, projectName, refType,diseaseFlag,bonferroniFlag)

if ~exist('diseaseFlag','var')  diseaseFlag=0;  end
if ~exist('bonferroniFlag','var')           bonferroniFlag = 0;         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fooofFlag = 1;
knee_flag = 0;   % for evaluating exponent for subject average psd
%freq_range = [4,120];
powerType = 'BL';
%electrodeGroupNum = 4; %group list for plotting data

if fooofFlag
    if knee_flag
        str_condition = 'withKnee';
    else
        str_condition = 'withoutKnee';
    end
else
    str_condition = 'withoutFooof';
end

if ~exist('capType','var')  capType = 'Acticap64';  end
showSEMFlag=1;
medianFlag = 1;
OptimR_SQ = 0;
OptimExponent = 0.01;
smoothSigma = 5;
parameter = 'r_SQ';
cLimParam = [0 4]; % limit of plots
nanValue =0; %the value of parameter
groupNum = 0; %electrode group num to be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%variable parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tapers = 1;%1/2/3[1 2 3];
freqRangeWidths = [100 76 200];
freqRangeList = {[64 140],[230 430]};
freqList = [{'LFR (64-140 Hz)'} {'HFR (230-430 Hz)'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the High priority electrodes
%need to be independant

if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
 [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    groupNameList{length(groupNameList)+1} = 'highPriority';
end


%%%%%%%%%%%%%%%% for different groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numGroups = length(subjectNameList);

y=[];
err= [];
if length (folderSourceString)<=3
    totDataType = length(folderSourceString);
else
    totDataType = 1;
    folderSourceStringNew{1} = folderSourceString;
    clear folderSourceString;
    folderSourceString = folderSourceStringNew;
end

    dataTypeList{1} = 'Eyes Open'; dataTypeList{2} = 'Eyes Closed';

%if diseaseFlag
    for jDat = 1:totDataType

    if jDat == 2
        if diseaseFlag
         %for removing case subjects
            subjectNameList{1}([1 3 4]) = [];
            subjectNameList{2}([1 3 4]) = [];
        end
    end
    
    for i = 1:numGroups
        %[ageListnew{i}, genderListnew{i}] = getDemographicDetails(projectName, subjectNameList{i});
    %    [ ~,paramData{i}, powerElectrodeGroup{i},paramElectrodeGroupData{i},stdElectrodeGroupData{i},SubjectElectrodeGroupData{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freq_range, powerType, OptimR_SQ, OptimExponent,medianFlag);
        %for continuous freq range
         [freqVals{jDat},centreFreq,PSD{jDat}{i},paramData{jDat}{i}, powerElectrodeGroup{jDat}{i},paramElectrodeGroupData{jDat}{i},paramSubjectElectrodeGroupData{jDat}{i},logPSD{jDat}{i}, logPowerElectrodeGroup{jDat}{i},subjectNum{jDat}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(1),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
         %for Low freqency range (4-120)
        % [~,centreFreqSingle{1},~,paramDataSing{jDat}{1}{i}, ~,paramElectrodeGroupDataSing{jDat}{1}{i},paramSubjectElectrodeGroupDataSing{jDat}{1}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(2),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
         %for high Freq range (170-470)
         %[~,centreFreqSingle{2},~,paramDataSing{jDat}{2}{i}, ~,paramElectrodeGroupDataSing{jDat}{2}{i},paramSubjectElectrodeGroupDataSing{jDat}{2}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(3),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
%          MalePos{i} = strcmp('M',genderListnew{i});
%         subMale{i} = subjectNameList{i}(MalePos{i});
%          FemalePos{i} = strcmp('F',genderListnew{i});
%         subFemale{i} = subjectNameList{i}(FemalePos{i});
    end
    end
% else
%     load(fullfile(pwd,'result_plots','VarWithFreqRange','02-02-23','fig1Data.mat'))
% end

displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
figure;
colormap jet;
    colorNames = [1 0 0; 0 0 1; 1 0 1;0.49 0.18 0.56; 0.5 0.5 0.5; 0 1 1 ; 0 1 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%hot(8); %colorNames([1:3,end-2:end],:) = [];
    displaySettings.colorNames = colorNames;
    color1 = [1 0.75 0.13];    color2=[1 0.3 0]; %add colors for line

%% figures psd and exponent
 axTopo = [];    axParam=[]; axPSD = []; yLimPSD = [];  yLimParam = [];

 %hPSD = getPlotHandles(totDataType,1,[.1 .1 0.4 0.8],0.05, 0.05);   linkaxes(hPSD);

for jDat = 1:totDataType % data types-baseline and eyes closed
    centreFreqNew = centreFreq;
     clear dataPSD dataParam dataTopo
     for i = 1:numGroups
         if groupNum ==0
             if medianFlag
                dataPSD{i} = squeeze(nanmedian(logPSD{jDat}{i},2));
                dataParam{i} = squeeze(nanmedian(paramData{jDat}{i},2));
             else
                 dataPSD{i} = squeeze(nanmean(logPSD{jDat}{i},2));
                dataParam{i} = squeeze(nanmean(paramData{jDat}{i},2));
             end
         else
            dataPSD{i} = logPowerElectrodeGroup{jDat}{i}{groupNum};
            dataParam{i} = paramSubjectElectrodeGroupData{jDat}{i}{groupNum};
         end
     
        % dataParam{i}(:,2:2:length(centreFreq)) = [];
     end
%     dataDiffLog{jDat} = 10*(nanmedian(logPowerElectrodeGroup{jDat}{2}{groupNum})-nanmedian(logPowerElectrodeGroup{jDat}{1}{groupNum})); %in db
%     dataDiffLinear{jDat} = 100*(nanmedian(powerElectrodeGroup{jDat}{2}{groupNum})-nanmedian(powerElectrodeGroup{jDat}{1}{groupNum}));
%     
    hPlot(1)  = subplot(3,2,jDat);%subplot('Position',[.1 .58-.45*(jDat-1) 0.25 0.35]); %psd
    hPlot(2) = subplot(3,2,jDat+2);%subplot('Position',[0.4 0.58-0.45*(jDat-1) 0.3 0.35]);%exponent for jDat1
    
   % hPlot(3) = subplot('Position',[0.45 0.33 0.45 0.15]); % exponent for jDat2
    %centreFreqNew(2:2:length(centreFreq)) = [];
    displayAndcompareData(hPlot(1),dataPSD,freqVals{jDat},displaySettings,'',0,medianFlag); xlim([6 800]); 
    %displayAndcompareData(hPlot(1),dataPSD,freqVals{jDat},displaySettings,'',1,medianFlag,'','','',bonferroniFlag); xlim([6 800]); 
    displayAndcompareData(hPlot(2),dataParam,centreFreqNew,displaySettings,[0 3],1,medianFlag,'','','',bonferroniFlag);  xlim([6 800]);
     
    
    yLimPSD = [yLimPSD hPlot(1).YLim];
    yLimParam = [yLimParam hPlot(2).YLim];
    set(hPlot(1),'Xscale','log','LineWidth',1,'FontSize',10,'FontWeight','bold');  
    set(hPlot(2),'Xscale','log','LineWidth',1,'FontSize',10,'FontWeight','bold'); 
    
   
         ylabel(hPlot(1),'log_{10} (Power (\muV)^2)','FontWeight','bold','FontSize',12);
        ylabel(hPlot(2),'Slope','FontWeight','bold','FontSize',12);
        %xlabel(hPlot(1),'Frequency (Hz)');
      t = title(hPlot(1),dataTypeList{jDat},'FontSize',14,'Fontweight','bold');
     % xlabel(hPlot(2),'Frequency (Hz)','FontWeight','bold','FontSize',12);
      
%                titPos = t.Position;
   % t.Position = [0.5 0.1 t.Position(3)];
%                t.Position = [-0.8 titPos(2)-0.5 titPos(3)];
    
     

    
     axPSD = [axPSD hPlot(1)];
     axParam = [axParam hPlot(2)];

 
  
end

linkaxes(axPSD);
linkaxes(axParam);

 yRange1 = [min(yLimPSD) max(yLimPSD)];
 yRange2  = [min(yLimParam) max(yLimParam)];

 for jDat = 1:totDataType
%     makeBox(axPSD(jDat),freqRangeList{1},yRange1,color1 ,2,'--','V');
%     makeBox(axPSD(jDat),freqRangeList{2},yRange1,color2 ,2,'--','V');
      makeShadedRegion(axPSD(jDat),freqRangeList{1},yRange1,color1);
     makeShadedRegion(axPSD(jDat),freqRangeList{2},yRange1,color2);

     legend(axPSD(jDat),'','','',[strList{1} '(' num2str(subjectNum{jDat}{1}) ')'],'',[strList{2} '(' num2str(subjectNum{jDat}{2}) ')'],'FontSize',10,'Location','southwest');

%     makeBox(axParam(jDat),freqRangeList{1},yRange2,color1 ,2,'--','V');
%     makeBox(axParam(jDat),freqRangeList{2},yRange2,color2 ,2,'--','V');
      makeShadedRegion(axParam(jDat),freqRangeList{1},yRange2,color1);
     makeShadedRegion(axParam(jDat),freqRangeList{2},yRange2,color2);
 end

 %%%%%%%%%%%%%%difference plots%%%%%%%%%%%%
hPlot(3) = subplot(3,2,5); %difference in log scale
    hPlot(4) = subplot(3,2,6); %difference in linear scale
   
 for jDat=totDataType:-1:1
    h1(jDat) =  plot(hPlot(3),freqVals{jDat},smooth(dataDiffLog{jDat},smoothSigma),'color',colorNames(jDat+3,:),'LineWidth',1.5); hold(hPlot(3),'on');
    h2(jDat) = plot(hPlot(4),freqVals{jDat},smooth(dataDiffLinear{jDat},smoothSigma),'color',colorNames(jDat+3,:),'LineWidth',1.5);  hold(hPlot(4),'on');
 end
plot(hPlot(3),freqVals{jDat},0.*freqVals{jDat}.^0,'k','LineWidth',1);
    plot(hPlot(4),freqVals{jDat},0.*freqVals{jDat}.^0,'k','LineWidth',1); 
set(hPlot(3),'Xscale','log','LineWidth',1,'FontSize',10,'FontWeight','bold','TickDir','out','TickLength',displaySettings.tickLengthMedium,'Box','off','XLim',[6 800]);
    set(hPlot(4),'Xscale','log','LineWidth',1,'FontSize',10,'FontWeight','bold','TickDir','out','TickLength',displaySettings.tickLengthMedium,'Box','off','XLim',[6 800],'YLim',[-0.5 0.5]);
     for k = 3:4
%      makeBox(hPlot(k),freqRangeList{1},hPlot(k).YLim,color1 ,2,'--','V');
%     makeBox(hPlot(k),freqRangeList{2},hPlot(k).YLim,color2 ,2,'--','V');
        makeShadedRegion(hPlot(k),freqRangeList{1},hPlot(k).YLim,color1);
        makeShadedRegion(hPlot(k),freqRangeList{2},hPlot(k).YLim,color2);
     end
     xlabel(hPlot(3),'Frequency (Hz)','FontWeight','bold','FontSize',12);
     xlabel(hPlot(4),'Frequency (Hz)','FontWeight','bold','FontSize',12);
     ylabel(hPlot(3),'\Delta Power (dB)','FontWeight','bold','FontSize',12);
     ylabel(hPlot(4),'\Delta Power (10^{-2} * (\muV)^2)','FontWeight','bold','FontSize',12);
      l = legend(h2([1 2]),dataTypeList{1},dataTypeList{2});
      l.Position(1) = 0.85;
     
%linkaxes(axTopo);
%annotation('textbox',[.2, .87 .1 .1],'String','PSD','FontSize',10,'Fontweight','bold','EdgeColor','none');

annotation('textbox',[.9, .75 .1 .05],'String',freqList{1},'color',color1,'FontSize',10,'fontWeight','bold','EdgeColor','none');
annotation('textbox',[.9, .72 .1 .05],'String',freqList{2},'color',color2,'FontSize',10,'fontWeight','bold','EdgeColor','none');
if bonferroniFlag == 1
    str2 = [{'KW Test'},{'p<0.01'}];
    annotation('textbox',[.91, .55 .07 .1],'String',str2{1},'FontSize',10,'EdgeColor','none');
    annotation('textbox',[.91, .545 .07 .05],'String','p<0.05 (BF)','color','g','FontSize',10,'EdgeColor','none');
elseif bonferroniFlag ==2
    str2 = [{'KW Test'},{'p<0.05 (CC)'}];
    annotation('textbox',[.91, .55 .07 .1],'String',str2,'FontSize',10,'EdgeColor','none');
    annotation('textbox',[.91, .545 .07 .05],'String','p<0.01 (CC)','color','g','FontSize',10,'EdgeColor','none');
else
    str2 = [{'KW Test'},{'p<0.05'}];
    annotation('textbox',[.91, .55 .07 .1],'String',str2,'FontSize',10,'EdgeColor','none');
    annotation('textbox',[.91, .545 .07 .05],'String','p<0.01','color','g','FontSize',10,'EdgeColor','none');
end

a = annotation('textbox',[0.0407 0.9286 0.0228 0.0579],'String','A','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[0.4907 0.9286 0.0228 0.0579],'String','B','FontSize',16,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[0.0407 0.3286 0.0228 0.0579],'String','C','FontSize',16,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[0.4907 0.3286 0.0228 0.0579],'String','D','FontSize',16,'EdgeColor','none','FontWeight','bold');
%if groupNum==0 sgtitle('ALL Electrodes');    else    sgtitle(groupNameList{groupNum},'FontWeight','bold');   end

end

function makeShadedRegion(hPlot,x,y,color,faceAlphaVal)
    if ~exist('faceAlphaVal','var') faceAlphaVal=0.3;   end
    xx = [x(1) x(1) x(2) x(2)];
    yy = [y(1)  y(2) y(2) y(1)];
    patch(hPlot,xx,yy,color,'LineStyle','none','FaceAlpha',faceAlphaVal);
    h = get(hPlot,'Children');
    set(hPlot,'Children',circshift(h,length(h)-1));
end
