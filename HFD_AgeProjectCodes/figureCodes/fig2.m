%% Figure 2 of paper

folderName = fullfile(folderSourceString{1},projectName,'higFracDim');
stFlag = 1;
medianFlag=1;
diseaseFlag = 0;
dimToMean = 2;
clear dim0 dimS0 dimDiff parameter
for i=1:length(subjectNameList)
    numSubjects(i) = 0;
    for j=1:length(subjectNameList{i})
        fileName = fullfile(folderName,[subjectNameList{i}{j} '_' refType '.mat']);
        if ~isfile(fileName)
            disp(['fileName for subject ' subjectNameList{i}{j} ' does not exist']);
            for k=1:8
                dim0{k}{i}(j,1:64) = nan;
                dimS0{k}{i}(j,1:64) = nan;
            end

        else
            xj = load(fileName);
            numSubjects(i) = numSubjects(i)+1;

            for k=1:length(xj.fracDim.higBL2)
                dim0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higBL2{k}(:,:,1),dimToMean));
                if stFlag
                    dimS0{k}{i}(j,:) = squeeze(nanmean(xj.fracDim.higSt2{k}(:,:,1),dimToMean));
                end
            end


        end
    end
end

if stFlag
    for i=1:2
        for k=1:length(xj.fracDim.higBL2)
            dimDiff{k}{i} = dimS0{k}{i}-dim0{k}{i};
        end
    end
end

parameter{1} = dim0{8};
parameter{2} = dimS0{8};
parameter{3} = dimDiff{8};


%% plotting
if strcmp(refType,'bipolar')
    electrodeGroupList{1} = [94, 93, 101, 96, 97, 102, 107, 111, 112];
    groupNameList{1} = 'highPriority';
else
    [~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(64,'EEG');
    %         electrodeGroupList{length(electrodeGroupList)+1} = highPriorityElectrodeNums;
    %         groupNameList{length(groupNameList)+1} = 'highPriority';
end
for iNum = 1:length(groupNameList)
    if isempty(find(groupNameList{iNum}=='-',1))
        groupNameShort{iNum} = groupNameList{iNum}(1);
    else
        groupNameShort{iNum} = [groupNameList{iNum}([1 (find(groupNameList{iNum}=='-'))+1])];
    end
end
numGroup = length(subjectNameList);

capType = 'actiCap64';
%% plotting the barplots for data


totCondns = length(parameter);
%medianFlag=1;

for iCond = 1:totCondns
    dataY{iCond} = [];
    for iElecGroup = 1:length(groupNameList)
        for iSub = 1:numGroup
            if medianFlag
                dataViol{iCond}{iElecGroup}{iSub} = nanmedian(parameter{iCond}{iSub}(:,electrodeGroupList{iElecGroup}),2);
                dataBar{iCond}(iElecGroup,iSub) = nanmedian(nanmedian(parameter{iCond}{iSub}(:,electrodeGroupList{iElecGroup}),2));
                dataTopo{iCond}{iSub} = nanmedian(parameter{iCond}{iSub},1);
            else
                dataViol{iCond}{iElecGroup}{iSub} = nanmean(parameter{iCond}{iSub}(:,electrodeGroupList{iElecGroup}),2);
                dataTopo{iCond}{iSub} = nanmean(parameter{iCond}{iSub},1);
            end
            dataY{iCond} = [dataY{iCond} dataViol{iCond}{iElecGroup}{iSub}'];
        end
    end
end

%for Topoplots
if strcmp(refType,'unipolar')
    cL = load([capType '.mat']);
    chanlocs = cL.chanlocs;
else
    cL = load(['bipolarChanlocs' capType '.mat']);
    chanlocs = cL.eloc;
end
nanValue = 0;


%% Plots
hf = figure;
hf.Position =  [142  166  1075  627];
%  sgtitle([num2str(freqRange{1}(1)) '-' num2str(freqRange{1}(2))]);
%hPlotViol = getPlotHandles(totCondns,1,[0.1,0.07,0.4,0.9],0.05,0.05);
hPlot = getPlotHandles(totCondns,3,[0.53,0.07,0.42,0.9],0.02,0.01);
colormap("jet")
colorNamesBl = [1 0 0; 0 0 1];

elecgroupTot = length(electrodeGroupList);
colorNames = copper(2*elecgroupTot);
colorNames(1:elecgroupTot,:) = colorNames(2:2:2*elecgroupTot,:);

markertype = [{'d'},{'v'},{'s'},{'+'},{'o'}];

diffParam = [0.01 0.01 0.001];%0.0005; % used to set ylims
tempx = groupNameShort;
tempRep = repelem(tempx,numSubjects(1)+numSubjects(2));
posY = [0.65 0.33 0.007];
yLabels = [{'HFD'},{'HFD'}, {'HFD'}];
clear g

for iCond = 1:totCondns

    %% 1. bar plots
    plotNumBar = iCond;%4*(iCond-1)+1;%
    % violpin plots

    violdata.x = tempRep;
    violdata.y =dataY{iCond} ;
    z1 = repelem(strList,[numSubjects(1) numSubjects(2)])';
    z = repelem(z1,1,length(groupNameList));
    violdata.group = z(:);
    g(1,iCond)=gramm('x',violdata.x,'y',violdata.y,'color',violdata.group);
    g(1,iCond).stat_violin('fill','transparent','normalization','width');
    g(1,iCond).stat_boxplot('width',0.13);
    if iCond==3;    g.geom_hline('yintercept',0,'style','--k');         end
    % g(1,iCond).set_title(titleName{iCond},'FontSize',15);
    g(1,iCond).set_color_options('map',colorNamesBl);
    g(1,iCond).axe_property('FontWeight','bold','FontSize',13,'TickLength',[0.01 0.01]);
    g(1,iCond).set_layout_options('position',[0.04 posY(iCond) 0.49 0.28],'legend_position',[0 0 0 0]);
    g(1,iCond).set_names('x','','y',yLabels{iCond},'Color','');
    g(1,iCond).set_order_options('x',0,'color',0);
    if iCond~=totCondns
        g(1,iCond).axe_property('xTickLabels',[],'FontWeight','bold','FontSize',13);
    else
        g(1,iCond).axe_property('FontWeight','bold','FontSize',13);
    end

    %% 2.topoplots
    %2.1 individual plots

    clear caxis
    for iSub = 1:numGroup
        plotNumTopo = (totCondns)*(iSub-1)+iCond;%4*(iCond-1)+iSub+1; %
        axes(hPlot(plotNumTopo));
        dataTopo{iCond}{iSub}(isnan(dataTopo{iCond}{iSub})) = nanValue;
        topoplot_murty(dataTopo{iCond}{iSub},chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir','+X','emarkercolors', dataTopo{iCond}{iSub});
        % topoplot(dataTopo{iCond}{iSub},chanlocs,'emarkercolors', dataTopo{iCond}{iSub});
        %c = colorbar;   c.Location = 'southoutside';
        if iCond~=3
            caxis([1.08 1.12]);
        else
            caxis([-0.002 0.002])
        end

        if iCond ==1;    title(strList{iSub},'FontSize',14);       end

        if iSub==1
            if iCond==1
                c=colorbar;
                c.Location = 'southoutside';
                c.Position = [0.6038 0.3803 0.1421 0.0270];
                c.FontSize = 10; c.FontWeight = 'bold';

                topoplot_murty([],chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','plotchans',highPriorityElectrodeNums,'emarker',{'x','k',12,1});
            elseif iCond==3
                c=colorbar;
                c.Location = 'southoutside';
                c.Position = [0.6038 0.07 0.1421 0.0270];
                c.FontSize = 10; c.FontWeight = 'bold';
            end
        end
    end

    %2.2 diff topoplots
    plotNumDiff = 2*totCondns+iCond;%4*iCond;%
    axes(hPlot(plotNumDiff));
    dataDiff{iCond} = dataTopo{iCond}{2} - dataTopo{iCond}{1};
    topoplot_murty(dataDiff{iCond},chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir','+X','emarkercolors', dataDiff{iCond});

    hold on;

    sigElectrodes = findSignificantElectrodes(parameter{iCond},medianFlag,0);
    if length(find(isnan(sigElectrodes)))~=64
        topoplot_murty([],chanlocs,'electrodes','on','style','blank','drawaxis','off','nosedir','+X','plotchans',find(sigElectrodes==1),'emarker',{'.','k',12,1})
    end

    if iCond~=totCondns;  caxis([-.02 .02]); else   caxis([-0.002 0.002]);  end%caxis([-.002 .002]);  end

    if iCond==1
        title([strList{2} '-' strList{1}],'FontSize',14);
        c=colorbar;
        c.Location = 'southoutside';
        c.Position = [0.84  0.3803  0.1 0.0270];
        c.FontSize = 10; c.FontWeight = 'bold';
    elseif iCond==3
        c=colorbar;
        c.Location = 'southoutside';
        c.Position = [0.84 0.07 0.1 0.0270];
        c.FontSize = 10; c.FontWeight = 'bold';
    end

end
g.draw();

allAxesInFigure = findall(gcf,'type','axes');
gPosAxes = [6 4 2];
gTitleRot = [{'BL'}, {'ST'},{'ST-BL'}];
gTitleDiff = [0.062 0.062 0.004];
legendNum = [2 31];
% pvalue
for iCond = 1:totCondns
    h = allAxesInFigure(gPosAxes(iCond));
    tPosition = h.YLim(1)+gTitleDiff(iCond);
    pPosition = h.YLim(2)+diffParam(iCond)/2;
    for iElecGroup = 1:length(groupNameList)
        if medianFlag
            pD(iCond,iElecGroup)=ranksum(nanmedian(parameter{iCond}{1}(:,electrodeGroupList{iElecGroup}),2),nanmedian(parameter{iCond}{2}(:,electrodeGroupList{iElecGroup}),2));
        else
            %t-test
            [~,pD(iCond,iElecGroup)]=ttest2(nanmean(parameter{iCond}{1}(:,electrodeGroupList{iElecGroup}),2),nanmean(parameter{iCond}{2}(:,electrodeGroupList{iElecGroup}),2));
        end

        % displaying p values
        %         pPosition = yLim(2) - diffParam(iCond)/5;
        % if pD(iCond,iElecGroup)<0.05
        if pD(iCond,iElecGroup)<0.001
            text(h,iElecGroup-0.2,pPosition,num2str(pD(iCond,iElecGroup),'%.1e'),'FontSize',10,'FontWeight','bold');
        else
            text(h,iElecGroup-0.2,pPosition,num2str(pD(iCond,iElecGroup),'%.3f'),'FontSize',10,'FontWeight','bold');
        end
    end
    % end
    text(h,-1.2,tPosition,gTitleRot{iCond},'FontSize',16,'FontWeight','bold','Rotation',90);
    if iCond==totCondns
        text(h,-0.6824,-0.0098,'\Delta','FontSize',14,'FontWeight','bold','Rotation',90);
    end

    %add legend
    if iCond==1
        for i = setdiff(1:legendNum(2),legendNum);   legendData{i} = {''};    end
        for i=1:2;      legendData{legendNum(i)} = strList{i};  end
        [l1,l2] = legend(string(strList),'Location',[0.47 0.68 0.05 0.05],'Box','off','FontSize',12,'FontWeight','Bold');
        set(l2,'lineWidth',2)
        l2(3).XData(1) = 0.4;
        l2(5).XData(1) = 0.4;
    end

end
a = annotation('textbox',[.0125 0.9107 0.0228 0.0579],'String','(A)','FontSize',18,'EdgeColor','none','FontWeight','bold');
b = annotation('textbox',[.0125 0.5951 0.0228 0.0579],'String','(B)','FontSize',18,'EdgeColor','none','FontWeight','bold');
c = annotation('textbox',[.0125 0.2858 0.0228 0.0579],'String','(C)','FontSize',18,'EdgeColor','none','FontWeight','bold');


function [sigElecNum] = findSignificantElectrodes(data,useMedianFlag,fdrFlag)
if ~exist('fdrFlag','var');  fdrFlag=1;   end

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
