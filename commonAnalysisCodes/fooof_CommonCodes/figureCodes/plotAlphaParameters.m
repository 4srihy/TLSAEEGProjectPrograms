%%plot alpha parameters

function plotAlphaParameters(subjectNameList,strList,folderSourceString, projectName, refType,medianFlag,diseaseFlag)
if ~exist('medianFlag','var')   medianFlag = 1; end
if ~exist('diseaseFlag','var')   diseaseFlag = 0; end

parameters = [{'CF'}, {'PH'}, {'BW'}, {'AP'}, {'APOld'}];
parameterNames = [{'Centre Frequency'}, {'Peak Height'}, {'Bandwidth'}, {'Alpha Power above aperiodic activity'}, {'log_{10}(Alpha Power) full'}];

[~,~,~,electrodeGroupList] = electrodePositionOnGrid(64,'EEG');
numGroups = length(subjectNameList);
%data
for i=1:numGroups
    [g{i},subjectNum(i)] = getMeanGroupFooofDataVFreq(subjectNameList{i},folderSourceString, projectName, refType, 100, 'BL', 0, 0,0,1,'withoutKnee',1);
    psdOcc{i}=squeeze(log10(nanmedian(g{i}.SpecPower(:,electrodeGroupList{1},:),2))); 
end

%%figure 
figure;
%display settings
displaySettings.fontSizeLarge = 8; displaySettings.tickLengthMedium = [0.025 0];
colorNames = [1 0 1; 0 1 1 ; 0 1 0 ; 0 0 1 ; 1 0 0 ; 1 1 0 ; 0 1 0 ; 1 0 1];%hot(8); %colorNames([1:3,end-2:end],:) = [];
    displaySettings.colorNames = colorNames;

hPlot1 = subplot('Position',[0.33 0.65 0.3 0.3]);
hPlot2 = getPlotHandles(1,length(parameters),[0.05 0.1 0.9 0.4],0.05,0,0); %for parameters

%%%psd
 displayAndcompareData(hPlot1,psdOcc,g{2}.freq,displaySettings,'',0,medianFlag); xlim([4 800]);
    set(hPlot1,'Xscale','log','LineWidth',1,'FontSize',10,'FontWeight','bold');  
    ylabel(hPlot1,'log_{10} (Power (\muV)^2)','FontWeight','bold','FontSize',12);
    xlabel(hPlot1,'Frequency (Hz)','FontWeight','bold','FontSize',12);
     makeBox(hPlot1,[7 12],hPlot1.YLim,'k' ,1.5,'--','V');
    legend('',[strList{1} '(' num2str(subjectNum(1)) ')'],'',[strList{2} '(' num2str(subjectNum(2)) ')'],'FontSize',10);

 %%bar plots including swarm plots

  x = categorical(strList);
x = reordercats(x,strList);
axBar = []; axYlim = [];
 for j = 1:length(parameters)
     if j==5
         for i=1:numGroups
             g{i}.(parameters{j}) = log10(g{i}.(parameters{j})) ;
         end
     end

     axes(hPlot2(j));
        y=[];
        err = [];
        xAll = [];
        yAll = [];
        %axBar = [axBar hPlot(6)];
        for i=1:numGroups
            if medianFlag
                 stdParam{i}{j} = nanstd(bootstrp(10000,@nanmedian,g{i}.(parameters{j})));
                  y = [y nanmedian(g{i}.(parameters{j}))];
            else
               stdParam{i}{j} = nanstd(g{i}.(parameters{j})/sqrt(length(g{i}.(parameters{j}))));
                y = [y nanmean(g{i}.(parameters{j}))];
            end
           
            err = [err stdParam{i}{j}];
            yAll = [yAll g{i}.(parameters{j})];
            x1 = repmat(strList{i},length(g{i}.(parameters{j})),1)';
            if ~diseaseFlag     xAll = [xAll x1];  end
        end

         %significance level testing
         if medianFlag
            [pD(j),h(j),stats{j}]=ranksum(g{1}.(parameters{j}),g{2}.(parameters{j}));
        else
            %t-test
            [h(j),pD(j),xc(j,:,:),stats{j}]=ttest2(g{1}.(parameters{j}),g{2}.(parameters{j}));
         end


       

         b = bar(x,y,0.5,'FaceAlpha',0.8,'FaceColor',[0.4940 0.1840 0.5560],'LineWidth',1);
         set(hPlot2(j),'TickLength',[0.04, 0.02],'TickDir','out','LineWidth',1,'FontWeight','bold','Box','off');
         ylabel(parameterNames{j}); 
     
         hold on
        if diseaseFlag
            for i=1:numGroups
                %swarmchart(x(i),paramSubjectElectrodeGroupDataSing{jFreq}{i}{j}','color',[1 0.5 0],'MarkerFaceAlpha',0.3);
                plot(i*(g{i}.(parameters{j})).^0,g{i}.(parameters{j}),'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k');
            end

            for u= 1:length(subjectNameList{1})
                plot([1 2],[g{1}.(parameters{j})(u) g{2}.(parameters{j})(u)],'k','LineWidth',0.8);
            end

            %pPosition = 3.1;
        else

            xs = categorical(cellstr(xAll'),strList);
            swarmchart(xs,yAll,20,[1 0.5 0],'filled','MarkerFaceAlpha',0.5);

            er = errorbar(x,y,err,err,'LineWidth',1.5);%'o','Color','r','MarkerFaceColor','w','LineWidth',1.5);  %'#CB4779'  
            er.Color =[0 0 0];% [0.7500 0.3250 0.0980];%'b';%'#CB4779';%[0 0 0];                            
            er.LineStyle = 'none'; 

          %  pPosition = 3.8;
        end
     ax = gca;
     ax.YLim(1) =0;
     pPosition = ax.YLim(2);
        if pD(j)<0.001
         text(1,pPosition,['p = ' num2str(pD(j),'%.1e')]);
        else
           text(1,pPosition,['p = ' num2str(pD(j),'%.3f')]); 
        end
 end
        

        
       
end
