%figure 1 new

function plotFigureUnipolarBipolar_v2(subjectNameList,strList,folderSourceString, projectName)

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
OptimExponent = 0;
parameter = 'exponent';
cLimParam = [0 4]; % limit of plots
nanValue =0; %the value of parameter
groupNum = 1; %electrode group num to be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%variable parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tapers = 1;%1/2/3[1 2 3];
freqRangeWidths = [100 76 200];
freqList = [{'Low (64-140 Hz)'} {'High (230-430 Hz)'}];
refTypes = [{'unipolar'}, {'bipolar'}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the High priority electrodes
%need to be independant

% 

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

    dataTypeList{1} = 'Baseline'; dataTypeList{2} = 'Eyes Closed';

    figure;
colormap jet;

for iRef = 1:length(refTypes)
    refType = refTypes{iRef};

    if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
    else
     [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    %     electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    %     groupNameList{length(groupNameList)+1} = 'highPriority';
    end
clear centreFreqSingle paramDataSing paramElectrodeGroupDataSing paramSubjectElectrodeGroupDataSing
    for jDat = 1:totDataType
        
        for i = 1:numGroups
            %[ageListnew{i}, genderListnew{i}] = getDemographicDetails(projectName, subjectNameList{i});
        %    [ ~,paramData{i}, powerElectrodeGroup{i},paramElectrodeGroupData{i},stdElectrodeGroupData{i},SubjectElectrodeGroupData{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freq_range, powerType, OptimR_SQ, OptimExponent,medianFlag);
            %for continuous freq range
           %  [freqVals{jDat},centreFreq,PSD{jDat}{i},paramData{jDat}{i}, powerElectrodeGroup{jDat}{i},paramElectrodeGroupData{jDat}{i},paramSubjectElectrodeGroupData{jDat}{i},logPSD{jDat}{i}, logPowerElectrodeGroup{jDat}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(1),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
             %for Low freqency range (4-120)
             [~,centreFreqSingle{1},~,paramDataSing{jDat}{1}{i}, ~,paramElectrodeGroupDataSing{jDat}{1}{i},paramSubjectElectrodeGroupDataSing{jDat}{1}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(2),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
             %for high Freq range (170-470)
             [~,centreFreqSingle{2},~,paramDataSing{jDat}{2}{i}, ~,paramElectrodeGroupDataSing{jDat}{2}{i},paramSubjectElectrodeGroupDataSing{jDat}{2}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString{jDat}, projectName, refType, freqRangeWidths(3),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
    %          MalePos{i} = strcmp('M',genderListnew{i});
    %         subMale{i} = subjectNameList{i}(MalePos{i});
    %          FemalePos{i} = strcmp('F',genderListnew{i});
    %         subFemale{i} = subjectNameList{i}(FemalePos{i});
        end
    end


%% figures psd and exponent
 axTopo = [];    axParam=[]; axPSD = [];

 
for jDat = 1:totDataType % data types-baseline and eyes closed
  hPlotExp{iRef,jDat} = getPlotHandles(2,numGroups,[.1+0.45*(iRef-1) .55-(0.49*(jDat-1)) 0.25 0.4],0.02, 0.02);  % linkaxes(hPl);
 hPlotDiff{iRef,jDat} = getPlotHandles(2,1,[0.37+0.45*(iRef-1) .55-(0.49*(jDat-1)) 0.12 0.4],0.05, 0.02);

    %% topoplots
    if strcmp(refType,'unipolar')
        cL = load([capType '.mat']);
        chanlocs = cL.chanlocs;
    else
        cL = load(['bipolarChanlocs' capType '.mat']);
        chanlocs = cL.eloc;
    end

   
     for jFreq = 1:2
       for iNum=1:numGroups+1
           %hPlot = subplot(4,6,6*((2*(jDat-1)+jFreq)-1)+3*(iRef-1)+iNum);
          
           if iNum~=3
                axes(hPlotExp{iRef,jDat}(jFreq+2*(iNum-1)));
                if medianFlag
                dataTopo = squeeze(nanmedian(paramDataSing{jDat}{jFreq}{iNum}));%- squeeze(nanmedian(paramDataSing{jDat}{jFreq}{1}));
                else
                dataTopo = squeeze(nanmean(paramDataSing{jDat}{jFreq}{iNum}));%- squeeze(nanmean(paramDataSing{jDat}{jFreq}{1}));
                end

                topoplot(dataTopo,chanlocs,'electrodes','off','style','map','drawaxis','off','nosedir','+X','emarkercolors', dataTopo);
                 if iRef==1 
                     if jDat==1
                         caxis([0.5 2.5]);   
                     else
                         caxis([0.2 2]);
                     end
                 else
                     if jDat==1
                         caxis([0.5 2.3]);   
                     else
                         caxis([0.2 1.25]);
                     end
                 end
                if iNum==1 && iRef==1 
                    ht = text(-0.8,-0.5,freqList{jFreq},'FontWeight','bold');
                    set(ht,"Rotation",90);
                end
                if jDat==1 && jFreq ==1     title(strList{iNum}); end
           else
                axes(hPlotDiff{iRef,jDat}(jFreq))
             if medianFlag
                dataTopo = squeeze(nanmedian(paramDataSing{jDat}{jFreq}{2}))- squeeze(nanmedian(paramDataSing{jDat}{jFreq}{1}));
             else 
                dataTopo = squeeze(nanmean(paramDataSing{jDat}{jFreq}{2}))- squeeze(nanmean(paramDataSing{jDat}{jFreq}{1}));
             end
             topoplot_murty(dataTopo,chanlocs,'electrodes','off','style','map','drawaxis','off','nosedir','+X','emarkercolors', dataTopo); hold on;
           caxis([-0.5 0.5])
             sigElectrodes = findSignificantElectrodes(paramDataSing{jDat}{jFreq},medianFlag);
             if length(find(isnan(sigElectrodes)))~=64
              topoplot(sigElectrodes,chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','emarker',{'.','k',12,1});
             end
            if jDat==1 && jFreq ==1     title('Difference'); end
           end

           if (jFreq ==1) %&&  (jDat ==1)
              
               if iNum == 1
                    c1 = colorbar;
                   c1.Location = 'southoutside'; 
                   c1.Position = [0.185+0.45*(iRef-1) 0.52-0.49*(jDat-1) 0.08 0.02];
                    %if jDat == length(dataTypeList)    c1.Label.String = parameter;  end
               elseif iNum == 3
                    c1 = colorbar;
                    c1.Location = 'southoutside'; 
                    c1.Position = [0.39+0.45*(iRef-1) 0.52-0.49*(jDat-1) 0.08 0.02];
                   %if jDat == length(dataTypeList)  c1.Label.String = 'exponent difference';    end
               end
               c1.Label.FontWeight = 'bold';
               c1.Label.FontSize = 10;
               c1.FontSize = 8;
               c1.FontWeight='bold';
           end

         
        % hold on;
        
         %title(freqList{jFreq});
%         title([strList{i} 'Females (N=' num2str(length(subFemale{i})) ')']);
%         caxis(cLimParam);
%            if i==1  
%                t = title(hPlot(5),freqList{jFreq},"Rotation",90);
%                titPos = t.Position;
%                t.Position = [-0.8 titPos(2)-0.5 titPos(3)];
%            end
          % axTopo = [axTopo hTopo{jDat}(jFreq)];
          %colorbar;
      end
    end
end
a = annotation('textbox',[0.03 0.925 0.0228 0.0579],'String','(A)','FontSize',13,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[0.03 0.45 0.0228 0.0579],'String','(B)','FontSize',13,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[0.53 0.925 0.0228 0.0579],'String','(C)','FontSize',13,'EdgeColor','none','FontWeight','bold');
d = annotation('textbox',[0.53 0.45 0.0228 0.0579],'String','(D)','FontSize',13,'EdgeColor','none','FontWeight','bold');

end
%linkaxes(axTopo);
% annotation('textbox',[.01, .67 .1 .1],'String','Baseline','FontSize',13,'Fontweight','bold','EdgeColor','none');
% annotation('textbox',[.01, .22 .1 .1],'String','Eyes Closed','FontSize',13,'Fontweight','bold','EdgeColor','none');

% annotation('textbox',[.27, .9 .1 .1],'String','Unipolar','FontSize',13,'Fontweight','bold','EdgeColor','none');
% annotation('textbox',[.69, .9 .1 .1],'String','Bipolar','FontSize',13,'Fontweight','bold','EdgeColor','none');

annotation('textarrow',[0.03 0.08],[ 0.7 0.9],'String','Baseline', 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',18, 'color','k','FontWeight','bold', 'TextRotation',90);

annotation('textarrow',[0.03 0.08],[ 0.2 0.4],'String','Eyes Closed', 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',18, 'color','k','FontWeight','bold', 'TextRotation',90);
% c = colorbar;   c.Position =  [0.9 0.3 0.015 .4] ;
% c.Label.String = ['Exponent difference (' strList{2} '-'  strList{1} ')'];
% c.Label.FontWeight = 'bold';
% c.Label.FontSize = 10;
% annotation('textbox',[.2, .87 .1 .1],'String','PSD','FontSize',10,'Fontweight','bold','EdgeColor','none');
% 
% if groupNum==0 sgtitle('ALL Electrodes');    else    sgtitle(groupNameList{groupNum});   end

end
%% plotting electrode with significance in topoplot

function [sigElecNum] = findSignificantElectrodes(data,useMedianFlag)
    
    allData = [];   allIDs = [];
    for i=1:size(data,2)
         allData = cat(1,allData,data{i});
        allIDs = cat(1,allIDs,i+zeros(size(data{i},1),1));
    end
    clear p sigElecNum
    for j = 1:size(allData,2)
         if useMedianFlag
           p=kruskalwallis(allData(:,j),allIDs,'off');
         else
               [~,p]=ttest2(data{1}(:,j),data{2}(:,j)); % only tests 2 groups
         end

         if p<0.05
             sigElecNum(j) =1;
         else
             sigElecNum(j) = nan;
         end
    end
end