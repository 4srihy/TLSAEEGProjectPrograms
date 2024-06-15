%supplementary fig 1
%containing HFD, H and HFD+H
colorNamesBl = [1 0 1; 0 1 1; 0.5 0.8 0.7];
lineColor = [0.5 0.5 0.5];%[1 0.9 0.5];

%% to plot hfd and HI for all subjects
folderName = fullfile(folderSourceString{1},projectName,'higFracDim');
stFlag = 1;
medianFlag=1;
diseaseFlag = 0;
dimToMean = 2;
groupNum = 6;
clear dim0 dimS0 dimDiff parameter hurst0 hurstS0
for i=1:length(subjectNameList)
    if diseaseFlag
        if i==1;     medianFlag=1;    else;    medianFlag = 0; end
    end
    for j=1:length(subjectNameList{i})
        
        fileName = fullfile(folderName,[subjectNameList{i}{j} '_' refType '.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
            for k=1
                dim0{k}{i}(j,1:64) = nan;
                dimS0{k}{i}(j,1:64) = nan;
                hurst0{k}{i}(j,1:64) = nan;
                hurstS0{k}{i}(j,1:64) = nan;
            end
           
        else
            xj = load(fileName);

            for k=1:length(xj.fracDim.higBL2)
                dim0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higBL2{k}(:,:,1),dimToMean));
                dimS0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higSt2{k}(:,:,1),dimToMean));
                if k==1
                    hurst0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.hurst2,dimToMean));
                    hurstS0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.hurstSt2,dimToMean));
                    sum0{k}{i}(j,:) = dim0{k}{i}(j,:) + hurst0{k}{i}(j,:);
                end
                
            end
        end
    end
end

parameter{1}{1} = dim0{8};     parameter{1}{2} = dimS0{8};       parameterName{1} = 'HFD';
parameter{2}{1} = hurst0{1};   parameter{2}{2} = hurstS0{1};      parameterName{2} = 'HI';

%% plotting
if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
    [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    groupNameList{length(groupNameList)+1} = 'highPriority';
end

numGroup = length(subjectNameList);
capType = 'actiCap64';

totCondns = length(parameter);
iElecGroup = groupNum;

for iCond = 1:totCondns
    for iSub = 1:2
        if medianFlag
            dataSwarm{iCond}{iSub} = cat(1,nanmedian(parameter{iCond}{iSub}{1}(:,electrodeGroupList{iElecGroup}),2),nanmedian(parameter{iCond}{iSub}{2}(:,electrodeGroupList{iElecGroup}),2));
        else
            dataSwarm{iCond}{iSub} = cat(1,nanmean(parameter{iCond}{iSub}{1}(:,electrodeGroupList{iElecGroup}),2),nanmean(parameter{iCond}{iSub}{2}(:,electrodeGroupList{iElecGroup}),2));
        end
    end
    dataY{iCond} = [dataSwarm{iCond}{1}' dataSwarm{iCond}{2}'];
end
dataSwarm{3}{1} = dataSwarm{1}{1}+dataSwarm{2}{1};
dataSwarm{3}{2} = dataSwarm{1}{2}+dataSwarm{2}{2};
dataY{3} = [dataSwarm{3}{1}'  dataSwarm{3}{2}'];



%%  add violin plots  % legend off is not working

for iCond = 1:3
    numDataPoints(iCond) = size(dataY{iCond},2)/2;
    dataDiff{iCond} = dataY{iCond}(numDataPoints(iCond)+1:end)-(dataY{iCond}(1:numDataPoints(iCond)));
end


posx = [0.01, 0.327 , 0.69];
yLabels = [{'HFD'}, {'H'}, {'HFD+H'}];

tempx = [{'All Subjects'}, {'All Subjects'},{'All Subjects'}];

clear g
for iCond = 1:3

    %data diff
    xDiff = repelem(tempx,[numDataPoints(iCond)]);
        violdata{1}.x = xDiff(numDataPoints(iCond)+1:2*numDataPoints(iCond));
   
    violdata{1}.y = dataDiff{iCond};
    violdata{1}.group = violdata{1}.x;

   
    for s = 1 %diff plots
   
            g(s,iCond)=gramm('x',violdata{s}.x,'y',violdata{s}.y,'color',violdata{s}.group);
            g(s,iCond).stat_violin('fill','transparent','normalization','width');
            g(s,iCond).stat_boxplot('width',0.13);
            %g(1,iCond).set_title(titleName{iCond},'FontSize',15);
            if s==1
               
                g(s,iCond).set_color_options('map',colorNamesBl(3,:));
            else
                g(s,iCond).set_color_options('map',colorNamesBl);
            end
            g(s,iCond).axe_property('FontWeight','bold','FontSize',14);
           
            if s~=2
            g(s,iCond).set_layout_options('position',[posx(s,iCond) 0.012 0.29 0.45],'legend_position',[0 0 0 0]);
            else
                g(s,iCond).set_layout_options('position',[posx(s,iCond) 0.012 0.33 0.45],'legend_position',[0 0 0 0]);
            end
            g(s,iCond).set_names('x','','y',yLabels{iCond},'Color','');
            g(s,iCond).set_order_options('x',0,'color',0);
        
    end
end

g.draw();

%% adding violin plots including pairing lines
displaySettings.parametricTest = 0;
displaySettings.showYTicks=1;
displaySettings.showXTicks=1;
displaySettings.commonYLim = 1;
displaySettings.xPositionText =0.8;
displaySettings.textFontSize = 10;
displaySettings.yPositionLine=0.15;
displaySettings.xTickLabels = [{'BL'} {'ST'}];
hPlot = getPlotHandles(1,3,[0.07 0.59 0.9 0.3750],0.1,0.11,0);

for iCond = 1:3
    h1 = hPlot(iCond);
    displaySettings.plotAxes = h1;

   
     ylabel(h1,yLabels{iCond});
 ax=displayViolinPlot_Sri(dataSwarm{iCond},[{colorNamesBl(1,:)} {colorNamesBl(2,:)}],1,1,0,1,displaySettings);
 set(ax,'FontWeight','bold','FontSize',14);
  switch iCond
%         case 1
%              displaySettings.setYLim=[1.04 1.18];
%         case 2
%              displaySettings.setYLim=[0.8 0.95];   
        case 3
             ylim([1.95 2.015]);   
    end
end

allAxesInFigure = findall(gcf,'type','axes');
%p -value
diffParam = -[0.005 0.005 0.005];
deltaXPos = [0.1 -0.035 0.15];
deltaYPos = [-0.006 -0.0015 -0.0025];
for iCond = 1:3
    % p value Here only for equal sets
    numData = size(dataY{iCond},2)/2;
    if medianFlag
        pD(iCond)=signrank(dataY{iCond}(1:numData),dataY{iCond}(numData+1:end));
    else
        %t-test
        [~,pD(iCond)]=ttest(dataY{iCond}(1:numData),dataY{iCond}(numData+1:end)); %from 'ttest
    end

    h =allAxesInFigure(2*(6-iCond)-1);% allAxesInFigure(2*(6-iCond));
    plot(h,[h.XLim(1) h.XLim(2)-0.2],[0 0],'--k','LineWidth',1);
    pPosition = h.YLim(2)-diffParam(iCond)/5;

    if pD(iCond)<0.001
        t =   text(h,0.8,pPosition,num2str(pD(iCond),'%.1e'),'FontSize',12,'FontWeight','bold');
    else
        text(h,0.9,pPosition,num2str(pD(iCond),'%.3f'),'FontSize',12,'FontWeight','bold');
    end
    text(h,deltaXPos(iCond),deltaYPos(iCond),'\Delta','FontSize',14,'FontWeight','bold','Rotation',90);

end
a = annotation('textbox',[0 0.9525 0.0228 0.0579],'String','(A)','FontSize',16,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.3033 0.9525 0.0228 0.0579],'String','(B)','FontSize',16,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.6402 0.9525 0.0228 0.0579],'String','(C)','FontSize',16,'EdgeColor','none','FontWeight','bold');
