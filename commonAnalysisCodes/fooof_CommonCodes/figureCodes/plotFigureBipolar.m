%figure 1 new

function plotFigureBipolar(subjectNameList,strList,folderSourceString, projectName,refType)

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
parameter = 'exponent';
cLimParam = [0 4]; % limit of plots
nanValue =0; %the value of parameter
groupNum = 1; %electrode group num to be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%variable parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tapers = 1;%1/2/3[1 2 3];
freqRangeWidths = [100 76 200];%[40 116 300]%
freqList = [{'LFR (64-140 Hz)'} {'HFR (230-430 Hz)'}];
refTypes = [{'unipolar'}, {'bipolar'}];

%%%%%%%%%%%%%%%% for different groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numGroups = length(subjectNameList);

y=[];
err= [];
totDataType = 1;
figure;
colormap jet;%turbo;

iRef = find(strcmp(refTypes,refType));
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
                 %for Low freqency range (4-120)
             [~,centreFreqSingle{1},~,paramDataSing{jDat}{1}{i}, ~,paramElectrodeGroupDataSing{jDat}{1}{i},paramSubjectElectrodeGroupDataSing{jDat}{1}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidths(2),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
             %for high Freq range (170-470)
             [~,centreFreqSingle{2},~,paramDataSing{jDat}{2}{i}, ~,paramElectrodeGroupDataSing{jDat}{2}{i},paramSubjectElectrodeGroupDataSing{jDat}{2}{i}] = getGroupWiseElectrodeGroupData(subjectNameList{i},parameter, electrodeGroupList,folderSourceString, projectName, refType, freqRangeWidths(3),tapers, powerType, OptimR_SQ, OptimExponent,medianFlag,str_condition);
        end
    end


%% figures psd and exponent
 axTopo = [];    axParam=[]; axPSD = [];

 
for jDat = 1:totDataType % data types-baseline and eyes closed
  

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
          
           if iNum~=3
               hPlot = subplot(2,3,3*(jFreq-1)+iNum+1);
                if medianFlag
                dataTopo = squeeze(nanmedian(paramDataSing{jDat}{jFreq}{iNum}));%- squeeze(nanmedian(paramDataSing{jDat}{jFreq}{1}));
             else 
                dataTopo = squeeze(nanmean(paramDataSing{jDat}{jFreq}{iNum}));%- squeeze(nanmean(paramDataSing{jDat}{jFreq}{1}));
                end

                topoplot_murty(dataTopo,chanlocs,'electrodes','off','style','blank','drawaxis','off','nosedir','+X','emarkercolors', dataTopo); 
                caxis([0.5 2]); 

                
                   
                
                if jDat==1 && jFreq ==1     
                    t = title(strList{iNum},'FontSize',12); 
                    t.Position(2) = 0.58;
                end
           else
               hPlot = subplot(2,3,3*(jFreq-1)+1);
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
            if jDat==1 && jFreq ==1     
                t = title('Difference','FontSize',12);
                t.Position(2) = 0.58;
            end
             ht = text(-0.8,-0.5,freqList{jFreq},'FontWeight','bold','FontSize',12);
                    set(ht,"Rotation",90);
           end

           if (jFreq ==1) &&  (jDat ==1)
              
               if iNum == 1
                    c1 = colorbar;
                   c1.Location = 'southoutside'; 
                   c1.Position = [0.45 0.08 0.4 0.03];
                  if strcmp(parameter,'exponent') 
                      c1.Label.String = 'Slope';
                  else
                      c1.Label.String = parameter;
                  end
               elseif iNum == 3
                    c1 = colorbar;
                    c1.Location = 'southoutside'; 
                    c1.Position = [0.15 0.08 0.2 0.03];
                   c1.Label.String = 'Slope Difference';
               end
               c1.FontSize = 8;
               c1.FontWeight = 'bold';
               c1.Label.FontWeight = 'bold';
               c1.Label.FontSize = 8;
           end
      end
    end
end
end

%% plotting electrode with significance in topoplot

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
             p(k) = p(k)*length(p)/(find(pAsc==p(k)));
         end
     end 
             sigElecNum(p<0.05) = 1;
end