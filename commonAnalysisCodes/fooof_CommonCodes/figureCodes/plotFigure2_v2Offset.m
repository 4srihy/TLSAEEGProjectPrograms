%holds for unipolar or average ref type
%v2 - 12-10-22
function plotFigure2_v2Offset(subjectNameList,strList,folderSourceString, projectName, refType,bonferroniFlag,diseaseFlag,subjectRemovalFlag)

if ~exist('diseaseFlag','var')  diseaseFlag=0;  end
if ~exist('subjectRemovalFlag','var')  subjectRemovalFlag=0;  end

if diseaseFlag
    if subjectRemovalFlag %for removing case subjects for eyes closed
        subjectNameList{1}([1 3 4]) = [];
        subjectNameList{2}([1 3 4]) = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fooofFlag = 1;
knee_flag = 0;   % for evaluating exponent for subject average psd
%freq_range = [4,120];
powerType = 'BL';
electrodeGroupNum = 4; %group list for plotting data

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
parameter = 'offset';
cLimParam = [0 4]; % limit of plots
nanValue =0; %the value of parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%variable parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tapers = 1;%1/2/3[1 2 3];
freqRangeWidths = [100 76 200];
freqRangeList = {[64 130],[230 430]};
freqList = [{'LFR (64-140 Hz)'} {'HFR (230-430 Hz)'}];%[{'Low (4-120 Hz)'} {'High (170-470 Hz)'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get the High priority electrodes
%need to be independant

if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
 [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
%     electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
%     groupNameList{length(groupNameList)+1} = 'highPriority';
end


%%%%%%%%%%%%%%%% for different groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numGroups = length(subjectNameList);

y=[];
err= [];
for i = 1:numGroups
    %[ageListnew{i}, genderListnew{i}] = getDemographicDetails(projectName, subjectNameList{i});
%    [ ~,paramData{i}, powerElectrodeGroup{i},paramElectrodeGroupData{i},stdElectrodeGroupData{i},SubjectElectrodeGroupData{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freq_range, powerType, OptimR_SQ, OptimExponent,medianFlag);
    %for continuous freq range
     [freqVals,centreFreq,PSD{i},paramData{i}, powerElectrodeGroup{i},paramElectrodeGroupData{i},paramSubjectElectrodeGroupData{i},logPSD{i}, logPowerElectrodeGroup{i},subjectNum{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidths(1),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
     %for Low freqency range (4-120)
     [~,centreFreqSingle{1},~,paramDataSing{1}{i}, ~,paramElectrodeGroupDataSing{1}{i},paramSubjectElectrodeGroupDataSing{1}{i},subjectNumSing{1}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidths(2),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
     %for high Freq range (170-470)
     [~,centreFreqSingle{2},~,paramDataSing{2}{i}, ~,paramElectrodeGroupDataSing{2}{i},paramSubjectElectrodeGroupDataSing{2}{i},subjectNumSing{2}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidths(3),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
%      MalePos{i} = strcmp('M',genderListnew{i});
%     subMale{i} = subjectNameList{i}(MalePos{i});
%      FemalePos{i} = strcmp('F',genderListnew{i});
%     subFemale{i} = subjectNameList{i}(FemalePos{i});
end

%display settings
displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
figure;
mapColor = 'jet';
colormap(mapColor);
colorNames1 = jet(length(electrodeGroupList)); % for all electrode groups
    colorNames2 = [1 0 0; 0 0 1; 1 0 1; 0 1 1 ; 0 1 0 ; 0 0 1 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%general%hot(8); %colorNames([1:3,end-2:end],:) = [];
    %displaySettings.colorNames = colorNames;
color1 = [1 0.75 0.13];    color2=[1 0.3 0]; lineWidth = 1.25;%add colors for line
    xlimLow = 6;    xlimHigh = 800; %for PSD and Exponent
    if strcmp(parameter,'exponent')
        yLimsSet = [0 3];
        parameter = 'Slope';  %for plotting only
    elseif strcmp(parameter,'offset')
        yLimsSet = [-2 4];
        parameter = 'Offset'; 
    else
        yLimsSet = '';
    end

   
%%%% plots%%%%%%%%%%%%%%%%%%%
rTot = 4; %total rows
cTot = 7; %total columns

axPSD = [];  axParam=[];
hTopoDiff = getPlotHandles(2,1,[0.02 0.1 0.09 0.375], 0 ,0.05);
%% 1. For different electrode groups for each group
    displaySettings.colorNames = colorNames1;
    for i=1:numGroups
        hPlot(1) = subplot(rTot,cTot,i); %psds
        hPlot(2) = subplot(rTot,cTot,cTot+i); %slopes
        displayAndcompareData(hPlot(1),logPowerElectrodeGroup{i},freqVals,displaySettings,'',0,medianFlag); xlim([xlimLow xlimHigh]);
        displayAndcompareData(hPlot(2),paramElectrodeGroupData{i},centreFreq,displaySettings,'',0,medianFlag,'','','',bonferroniFlag); xlim([xlimLow xlimHigh]);
        set(hPlot(1),'XTickLabel',[],'Xscale','log','TickLength',[0.04 0.02],'XTick',[10 100],'LineWidth',1,'FontWeight','bold');
        set(hPlot(2),'Xscale','log','TickLength',[0.04 0.02],'XTick',[10 100],'LineWidth',1,'FontWeight','bold');
        hPlot(2).Position = hPlot(1).Position;
     hPlot(2).Position(2) = 0.59;
        xlabel(hPlot(2),'Frequency (Hz)');
        title(hPlot(1),strList{i},'FontSize',12);
        if i==1
            ylabel(hPlot(1),'log_{10}(Power (\muV^2))');
            ylabel(hPlot(2),parameter);
            if strcmp(refType,'unipolar')
           
            end
        end
        axPSD = [axPSD hPlot(1)];
        axParam = [axParam hPlot(2)];
    end

 %% 2. for different age groups for each electrode group
  yLimPSD = [];  yLimParam = [];
displaySettings.colorNames = colorNames2;
 for j = 1:length(electrodeGroupList)
%      hPlot(3) = subplot(rTot,cTot,cTot*(j+1)+1); %psd
%      hPlot(4) = subplot(rTot,cTot,cTot*(j+1)+2); %slope

    hPlot(3) = subplot(rTot,cTot,j+2); %psd
    hPlot(4) = subplot(rTot,cTot,cTot+(j+2)); %slope
     clear dataPSD dataParam
     for i = 1:numGroups
     dataPSD{i} = logPowerElectrodeGroup{i}{j};
     dataParam{i} = paramSubjectElectrodeGroupData{i}{j};
     end

     displayAndcompareData(hPlot(3),dataPSD,freqVals,displaySettings,'',0,medianFlag); xlim([xlimLow xlimHigh]);
     displayAndcompareData(hPlot(4),dataParam,centreFreq,displaySettings,yLimsSet,1,medianFlag,'','','',bonferroniFlag); xlim([xlimLow xlimHigh]);
      yLimPSD = [yLimPSD hPlot(3).YLim];
      yLimParam = [yLimParam hPlot(4).YLim];

     set(hPlot(3),'XTickLabel',[],'Xscale','log','TickLength',[0.04 0.02],'XTick',[10 100],'LineWidth',1,'FontWeight','bold'); 
     set(hPlot(4),'Xscale','log','TickLength',[0.04 0.02],'XTick',[10 100],'LineWidth',1,'FontWeight','bold');
     hPlot(4).Position = hPlot(3).Position;
     hPlot(4).Position(2) = 0.59;
      xlabel(hPlot(4),'Frequency (Hz)');
%      if j==length(groupNameList)  
%          l2 = legend(hPlot(3),'',strList{1},'',strList{2});
%          l2.Position = [0.91 0.85 0.07 0.05];
%      end
     
      title(hPlot(3),groupNameList{j},'FontSize',12);
      axPSD = [axPSD hPlot(3)];
        axParam = [axParam hPlot(4)];
 end
linkaxes(axPSD);
linkaxes(axParam);

 yRange1 = [min(yLimPSD) max(yLimPSD)];
 yRange2  = [min(yLimParam) max(yLimParam)];

 for j = 1:length(electrodeGroupList)+2
%     makeBox(axPSD(j),freqRangeList{1},yRange1,color1 ,lineWidth,'--','V');
%     makeBox(axPSD(j),freqRangeList{2},yRange1,color2 ,lineWidth,'--','V');
% 
%     % legend(axPSD(j),'',[strList{1} '(N=' num2str(subjectNum{j}{1}) ')'],'',[strList{2} '(N=' num2str(subjectNum{j}{2}) ')'],'FontSize',10);
% 
%     makeBox(axParam(j),freqRangeList{1},yRange2,color1 ,lineWidth,'--','V');
%     makeBox(axParam(j),freqRangeList{2},yRange2,color2 ,lineWidth,'--','V');
     makeShadedRegion(axPSD(j),freqRangeList{1},yRange1,color1);
     makeShadedRegion(axPSD(j),freqRangeList{2},yRange1,color2);

      makeShadedRegion(axParam(j),freqRangeList{1},yRange2,color1);
     makeShadedRegion(axParam(j),freqRangeList{2},yRange2,color2);

 end
 l= legend(hPlot(2),'','','',groupNameList{1},'',groupNameList{2},'',groupNameList{3},'',groupNameList{4},'',groupNameList{5});
            l.Position = [0 .8 .1 .1];
  l2 = legend(hPlot(3),'','','',[strList{1} '(' num2str(subjectNum{1}) ')'],'',[strList{2} '(' num2str(subjectNum{2}) ')'],'FontSize',9);
  
  l2.Position = [0.91 0.85 0.07 0.05];

annotation('textbox',[.91, .75 .1 .05],'String',freqList{1},'color',color1,'FontSize',9,'fontWeight','bold','EdgeColor','none');
annotation('textbox',[.91, .72 .1 .05],'String',freqList{2},'color',color2,'FontSize',9,'fontWeight','bold','EdgeColor','none');

if ~bonferroniFlag
    str2 = [{'KW Test'},{'p<0.05'}];
    annotation('textbox',[.91, .59 .07 .1],'String',str2,'FontSize',9,'EdgeColor','none');
    annotation('textbox',[.91, .585 .07 .05],'String','p<0.01','color','g','FontSize',9,'EdgeColor','none');
elseif  bonferroniFlag==1
    str2 = [{'KW Test'},{'p<0.01'}];
    annotation('textbox',[.91, .59 .07 .1],'String',str2,'FontSize',9,'EdgeColor','none');
    annotation('textbox',[.91, .585 .07 .05],'String','p<0.05 (BF)','color','g','FontSize',9,'EdgeColor','none');
elseif bonferroniFlag==2
     str2 = [{'KW Test'},{'p<0.05 (CC)'}];
    annotation('textbox',[.91, .59 .07 .1],'String',str2,'FontSize',9,'EdgeColor','none');
    annotation('textbox',[.91, .585 .07 .05],'String','p<0.01 (CC)','color','g','FontSize',9,'EdgeColor','none');
end
 %% 3. for Topoplots
 if strcmp(refType,'unipolar')
    cL = load([capType '.mat']);
    chanlocs = cL.chanlocs;
else
    cL = load(['bipolarChanlocs' capType '.mat']);
    chanlocs = cL.eloc;
 end

 %%%3.0 Elec topoplot
 hTopo = subplot('Position',[0.017 0.59 0.0864 0.1577]);
 elecVals = zeros(1,length(chanlocs));
 for iElec=1:length(electrodeGroupList)
     elecVals(electrodeGroupList{iElec}) = iElec;
 end

  ced = topoplot_murty([],chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','emarkercolors', elecVals); 
  %ced.colormap = 'summer';
   caxis([0 6.5]);

 %%%3.1 diff topoplots
  for jFreq = 1:2
         if medianFlag
            dataTopo = squeeze(nanmedian(paramDataSing{jFreq}{2}))- squeeze(nanmedian(paramDataSing{jFreq}{1}));
         else 
            dataTopo = squeeze(nanmean(paramDataSing{jFreq}{2}))- squeeze(nanmean(paramDataSing{jFreq}{1}));
         end
         axes(hTopoDiff(jFreq));
        % dataTopo(find(isnan(dataTopo))) = nanValue;
      % topoplot_murty(dataTopo,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir','+X','emarkercolors', dataTopo); hold on;
        topoplot_murty(dataTopo,chanlocs,'electrodes','off','style','map','drawaxis','off','nosedir','+X','emarkercolors', dataTopo); hold on;
        caxis([-1.5 1.5]);
       if strcmp(parameter,'exponent') caxis([-0.5 0.5]);
       elseif strcmp(parameter,'offset') caxis([-1.5 1.5]); end
         sigElectrodes = findSignificantElectrodes(paramDataSing{jFreq},medianFlag);
         if length(find(isnan(sigElectrodes)))~=64
          topoplot(sigElectrodes,chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','emarker',{'.','k',12,1});
         end
        
         %title(freqList{jFreq});
%         title([strList{i} 'Females (N=' num2str(length(subFemale{i})) ')']);
       % caxis(cLimParam);
           
               t = title(hTopoDiff(jFreq),freqList{jFreq},"Rotation",90);
               titPos = t.Position;
               t.Position = [-0.5 titPos(2)-0.5 titPos(3)];
  end
c = colorbar;   c.Location = 'southoutside'; c.Position =  [0.02 0.07 0.09 .02]; c.FontSize = 8; c.FontWeight = 'bold';
c.Label.String = [parameter ' Difference'];
c.Label.FontSize = 9;  c.Label.FontWeight = 'bold';

%%%3.2 median topoplots
 axTopo = [];
 for i=1:numGroups
     for jFreq = 1:2
         if medianFlag
         dataTopo = squeeze(nanmedian(paramDataSing{jFreq}{i}));
         else 
             dataTopo = squeeze(nanmean(paramDataSing{jFreq}{i}));
         end
         hPlot(5) = subplot(rTot,cTot,cTot*(jFreq+1)+i); %topoplots
         dataTopo(find(isnan(dataTopo))) = nanValue;
         topoplot_murty(dataTopo,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir','+X','emarkercolors', dataTopo);
%         title([strList{i} 'Females (N=' num2str(length(subFemale{i})) ')']);
%         caxis(cLimParam);
%            if i==1  
%                t = title(hPlot(5),freqList{jFreq},"Rotation",90);
%                titPos = t.Position;
%                t.Position = [-0.8 titPos(2)-0.5 titPos(3)];
%            end
           caxis([-1.5 3]);
     end
 end
%linkaxes(axTopo);
c = colorbar;   c.Location = 'southoutside'; c.Position =  [0.14 0.07 0.16 .02] ; c.FontSize = 8; c.FontWeight = 'bold';
c.Label.String = parameter; c.Label.FontSize = 10; c.Label.FontWeight = 'bold';

 %% barplots for each electrode group
 
 x = categorical(strList);
x = reordercats(x,strList);
axBar = []; axYlim = [];
for jFreq  =1:2
%if ~diseaseFlag      
    axBar = []; axYlim = [];    
%end
 for j = 1:length(electrodeGroupList)
     
        hPlot(6) = subplot(rTot,cTot,cTot*(jFreq+1)+(j+2)); %psd
        hPlot(6).Position(2) = -0.23*(jFreq-1)+0.29;
        hPlot(6).Position(4) = 0.2;
        y=[];
        err = [];
        xAll = [];
        yAll = [];
        axBar = [axBar hPlot(6)];
        for i=1:numGroups
            if medianFlag
                 stdElectrodeGroupDataSing{jFreq}{i}{j} = nanstd(bootstrp(10000,@nanmedian,paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}));
            else
               stdElectrodeGroupDataSing{jFreq}{i}{j} = nanstd(paramSubjectElectrodeGroupDataSing{jFreq}{i}{j})/sqrt(length(paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}));
            end
            y = [y paramElectrodeGroupDataSing{jFreq}{i}{j}];
            err = [err stdElectrodeGroupDataSing{jFreq}{i}{j}];
            yAll = [yAll paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}'];
            x1 = repmat(strList{i},length(paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}),1)';
            if ~diseaseFlag     xAll = [xAll, x1];  end
        end

         %significance level testing
         if medianFlag
            [pD(j,jFreq),h(j,jFreq),stats{j,jFreq}]=ranksum(paramSubjectElectrodeGroupDataSing{jFreq}{1}{j},paramSubjectElectrodeGroupDataSing{jFreq}{2}{j});
        else
            %t-test
            [h(j,jFreq),pD(j,jFreq),xc(j,jFreq,:,:),stats{j,jFreq}]=ttest2(paramSubjectElectrodeGroupDataSing{jFreq}{1}{j},paramSubjectElectrodeGroupDataSing{jFreq}{2}{j});
         end


       

         b = bar(x,y,0.5,'FaceAlpha',0.8,'LineWidth',1);
         set(hPlot(6),'TickLength',[0.04, 0.02],'YTick',[0:1:4],'TickDir','out','LineWidth',1,'FontWeight','bold','Box','off');
        if j==1 ylabel(parameter);  end
       if jFreq~=2
        set(gca,'XTickLabel',[]);
       end
        
         hold on
        if diseaseFlag
            for i=1:numGroups
                %swarmchart(x(i),paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}','color',[1 0.5 0],'MarkerFaceAlpha',0.3);
                plot(i*(paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}).^0,paramSubjectElectrodeGroupDataSing{jFreq}{i}{j},'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');
            end

            for u= 1:length(subjectNameList{1})
                plot([1 2],[paramSubjectElectrodeGroupDataSing{jFreq}{1}{j}(u) paramSubjectElectrodeGroupDataSing{jFreq}{2}{j}(u)],'k','LineWidth',0.8);
            end

            %pPosition = 3.1;
            if jFreq == 1    pPosition = 1.875;  else  pPosition=4;  end
        else

            xs = categorical(cellstr(xAll'),strList);
            swarmchart(xs,yAll,20,[1 0.5 0],'filled','MarkerFaceAlpha',0.3);

            er = errorbar(x,y,err,err,'LineWidth',1.5);%'o','Color','r','MarkerFaceColor','w','LineWidth',1.5);  %'#CB4779'  
            er.Color =[0 0 0];% [0.7500 0.3250 0.0980];%'b';%'#CB4779';%[0 0 0];                            
            er.LineStyle = 'none'; 

            if jFreq == 1    pPosition = 6;  else  pPosition=7;  end
        end

       
        if pD(j,jFreq)<0.001
         text(1,pPosition,['p = ' num2str(pD(j,jFreq),'%.1e')]);
        else
           text(1,pPosition,['p = ' num2str(pD(j,jFreq),'%.3f')]); 
        end

        

        
        yPerr = y+err;
        yBar = get(hPlot(6), 'YTick');
        axYlim = [axYlim hPlot(6).YLim(2)];
        

%         if pD(j,jFreq)<0.005       
%             plot([1.3 1.5 1.7], [1 1 1]*(max(y)+0.25), '*k'); 
%             plot([1 1.15],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%              plot ([1.85 2],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%              plot([1 1], [yPerr(1)+0.05 max(y)+0.25],'k','LineWidth',1);
%              plot([2 2], [yPerr(2)+0.05 max(y)+0.25],'k','LineWidth',1);
%         elseif pD(j,jFreq)<0.01  
%             plot([1.45 1.55], [1 1]*(max(y)+0.25), '*k');
%               plot ([1 1.3],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%              plot([1.7 2],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%               plot([1 1], [yPerr(1)+0.05 max(y)+0.25],'k','LineWidth',1);
%              plot([2 2], [yPerr(2)+0.05 max(y)+0.25],'k','LineWidth',1);
%         %elseif pD(j,jFreq)<0.05    plot([1 2], [1 1]*max(yBar)*1.2, '-k',1.5,max(yBar)*1.2, '*k');
%         elseif pD(j,jFreq)<0.05    
%             plot(1.5,(max(y)+0.25), '*k');
%               plot ([1 1.35],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%              plot ([1.65 2],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%               plot([1 1], [yPerr(1)+0.05 max(y)+0.25],'k','LineWidth',1);
%              plot([2 2], [yPerr(2)+0.05 max(y)+0.25],'k','LineWidth',1);
%         else
%             text(1.35,(max(y)+0.26),'n.s.');
%              plot ([1 1.3],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%              plot ([1.7 2],[(max(y)+0.25) (max(y)+0.25)],'k','LineWidth',1);
%               plot([1 1], [yPerr(1)+0.05 max(y)+0.25],'k','LineWidth',1);
%              plot([2 2], [yPerr(2)+0.05 max(y)+0.25],'k','LineWidth',1);
%         end
        hold off
     end
     linkaxes(axBar);
    ylim([-2  max(axYlim)+0.2]);
 end
% linkaxes(axBar);
% ylim([0  max(axYlim)+0.4]);
%str = [{'*   p<0.05'},{'**  p<0.01'},{'*** p<0.005'}];
%annotation('textbox',[.91, .05 .07 .3],'String',str,'FontSize',9,'EdgeColor','none');
a = annotation('textbox',[.0126 0.9352 0.0228 0.0579],'String','A','FontSize',13,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.33 0.9352 0.0228 0.0579],'String','B','FontSize',13,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.0126 0.48 0.0228 0.0579],'String','C','FontSize',13,'EdgeColor','none','FontWeight','bold');
%d = annotation('textbox',[.098 0.48 0.0228 0.0579],'String','(D)','FontSize',13,'EdgeColor','none','FontWeight','bold');
e = annotation('textbox',[.33 0.48 0.0228 0.0579],'String','D','FontSize',13,'EdgeColor','none','FontWeight','bold');

end

function [sigElecNum] = findSignificantElectrodes(data,useMedianFlag,fdrFlag)
if ~exist('fdrFlag','var')  fdrFlag=1;   end
    
    allData = [];   allIDs = [];
    for i=1:size(data,2)
         allData = cat(1,allData,data{i});
        allIDs = cat(1,allIDs,i+zeros(size(data{i},1),1));
    end
    clear p sigElecNum
    sigElecNum = nan(1,size(allData,2));
    for j = 1:size(allData,2)
         if useMedianFlag
           p(j)=kruskalwallis(allData(:,j),allIDs,'off');
         else
               [~,p(j)]=ttest2(data{1}(:,j),data{2}(:,j)); % only tests 2 groups
         end
    end
     if fdrFlag
         pAsc = sort(p,'ascend');
         
         for k = 1:length(pAsc)
             p(k) = p(k)*length(p)/(find(pAsc==p(k),1));
         end
     end
             sigElecNum(p<0.05) = 1;
end

function makeShadedRegion(hPlot,x,y,color,faceAlphaVal)
    if ~exist('faceAlphaVal','var') faceAlphaVal=0.3;   end
    xx = [x(1) x(1) x(2) x(2)];
    yy = [y(1)  y(2) y(2) y(1)];
    patch(hPlot,xx,yy,color,'LineStyle','none','FaceAlpha',faceAlphaVal);
     h = get(hPlot,'Children');
    set(hPlot,'Children',circshift(h,length(h)-1));
end