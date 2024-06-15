%To get the fooof data for different freq ranges for each electrode



%%%%%%%%getting subject details%%%%%%%%%%%

function saveFooofParamsIndividualSubjectsForRandomTrials(folderSourceString,subjectName,projectName,refType,protocolType, powerType, knee_flag,fooofPlotFlag,rmsPlotFlag,noneElec,fooofFlag)

if ~exist('knee_flag','var')    knee_flag=1;            end
%if ~exist('freq_range','var')   freq_range = [5,70];    end
if ~exist('fooofPlotFlag','var')    fooofPlotFlag=0;    end
if ~exist('rmsPlotFlag','var')      rmsPlotFlag = 0;    end
if ~exist('powerType','var')    powerType = 'BL';       end
%if ~exist('specPowerFlag','var')    specPowerFlag=0;    end
if ~exist('noneElec','var')     noneElec = [];          end
if ~exist('fooofFlag','var')    fooofFlag = 1;          end

numTrialsForRand = 40;

foldername = fullfile(folderSourceString,projectName);
load(fullfile(foldername,protocolType,[subjectName '_' refType '_stRange_250_750.mat']));

totalTrials = size(blPowerVsFreqTopo,3);
numTrials = randi([1 totalTrials],1,numTrialsForRand);
if strcmp(powerType,'BL')
     SpecPower= nanmean(blPowerVsFreqTopo(:,:,numTrials),3);
elseif strcmp(powerType,'ST')
    SpecPower = nanmean(stPowerVsFreqTopo(:,:,numTrials),3);
else
    disp("Specify the power type");
end

freqRes = freqVals(2)-freqVals(1);
allFreqRangeWidths = [100 76  200];%[100 116 300];%[40] %[100]
if allFreqRangeWidths(1)==100
    centreFreq = 40:20:960;
else
    centreFreq =20:10:950;
end %40:20:960;%20:10:950;
%for g=[3,10,13]   SpecPower(g,:) = nan; end

%%%%%to give multiple frequence ranges together
minFreqVal = 4; maxFreqVal = 950;
fullFreqRange = [minFreqVal,maxFreqVal];   %overall freq_range



freqsToAvoidDeclared{1} = [0 0];
awayFromCentre = 4;
for freqAvoidcentre = 50:50:950
    freqsToAvoidDeclared{freqAvoidcentre/50 +1} = [freqAvoidcentre-awayFromCentre freqAvoidcentre+awayFromCentre];
end
    %freqsToAvoidDeclared = {[0 0] [46 54] [96 104] [146 154] [196 204] [246 254] [296 304] [346 354] [396 404] [446 454] [496 504] [546 554] [596 604] [646 654] [696 704] [746 754] [796 804] [846 854] [896 904] [946 954]}; % Hz
for iFreqRangeWidth =1:length(allFreqRangeWidths)
    
     freqRangeWidth = allFreqRangeWidths(iFreqRangeWidth);%12(\pm 6),20, 24,32
    %only for high-low category
    if freqRangeWidth==116
        centreFreq = 62;
    elseif freqRangeWidth==300
        centreFreq = 320;
    elseif freqRangeWidth == 76
        centreFreq = 102;
    elseif freqRangeWidth == 200
        centreFreq = 330;
    end

    disp(['processing data for freqRangeWidth: ' num2str(freqRangeWidth)]);
    g=[];
    if iFreqRangeWidth==1    g=noneElec;    end
    
   

    clear freq_range freqLength
    for jFreq = 1:length(centreFreq)
        freq_range{jFreq} = [centreFreq(jFreq)-freqRangeWidth/2,centreFreq(jFreq)+freqRangeWidth/2];
        if freq_range{jFreq}(1)<minFreqVal
            freq_range{jFreq}(1) = minFreqVal;
        end
    
        if freq_range{jFreq}(2)>maxFreqVal
            freq_range{jFreq}(2) = maxFreqVal;
        end
    end
        %%to remove freq to avoid 
        %note that if both the frequencies of a given peak is present, then it 
        %can be sustained. Only if One of the frequencies are present, then no
        %peaks can be generated and fitting pose a problem. hence to be removed
       
    %     freqLength = freq_range{j}(1):freqRes:freq_range{j}(2);
    %             indicesToRemove = [];
    %                for iFreqAvoid = 1:length(freqsToAvoid)
    %                    if xor(ismember(freqsToAvoid{iFreqAvoid}(1),freqLength),ismember(freqsToAvoid{iFreqAvoid}(2),freqLength))
    %                         indicesToRemove{iFreqAvoid} = cat(2,indicesToRemove,find(freqLength>=freqsToAvoid{iFreqAvoid}(1) & freqLength<=freqsToAvoid{iFreqAvoid}(2)));
    %                    end
    %                end
    %            freqLength(indicesToRemove)  = nan;
    %            freq_rangeNew{j} = [nanmin(freqLength), nanmax(freqLength)];
    
    
    settings = struct();
    if strcmp(protocolType,'eyes_open')|| strcmp(protocolType,'eyes_closed')
         settings.peak_width_limits = [1,8];
    else
         settings.peak_width_limits = [4,8];
    end
    
    settings.max_n_peaks = 5;
    settings.min_peak_height = 0.2;
    if knee_flag      settings.aperiodic_mode = 'knee';      end
    %settings.peak_threshold = 1;
    fooof_results=[];
    exponent = zeros(size(SpecPower,1),length(freq_range)); 
    offset=zeros(size(SpecPower,1),length(freq_range)); 
    knee=zeros(size(SpecPower,1),length(freq_range)); 
    r_SQ=zeros(size(SpecPower,1),length(freq_range)); 
    icheck=[];
    icount=1;
    
     for ind =  1:size(SpecPower,1) %[24, 26, 29, 30, 31, 57, 58, 61, 62, 63]%   %index of the cases/control
        
         fooof_resultSingElec =[];
         if isnan(SpecPower(ind,1)) || ismember(ind,g)
             exponent(ind,:) = nan; offset(ind,:)=nan; knee(ind,:)=nan; r_SQ(ind,:)=nan;
         else
            freqsToAvoid = [];
            fooof_indFull = fooof(freqVals,SpecPower(ind,:)',fullFreqRange,settings,true);
            peakParams = fooof_indFull.peak_params;
%             disp(ind);
%             disp(peakParams);
            for iPeak = 1:size(peakParams,1)
                freqsToAvoid{iPeak} = [peakParams(iPeak,1)-peakParams(iPeak,3)/2 peakParams(iPeak,1)+peakParams(iPeak,3)/2];
            end
                freqsToAvoid = [freqsToAvoid freqsToAvoidDeclared];

            if fooofFlag   % computes using fooof    
             for j=1:length(freq_range)
                 freqLength = freq_range{j}(1):freqRes:freq_range{j}(2);
                indicesToRemove = [];
                   for iFreqAvoid = 1:length(freqsToAvoid)
                       if xor(length(find(freqLength<freqsToAvoid{iFreqAvoid}(1))),length(find(freqLength>freqsToAvoid{iFreqAvoid}(2))))
                         indicesToRemove = [indicesToRemove find(freqLength>=freqsToAvoid{iFreqAvoid}(1) & freqLength<=freqsToAvoid{iFreqAvoid}(2))];
                       end
                   end
               freqLength(indicesToRemove)  = nan;
               freq_rangeNew{j} = [nanmin(freqLength), nanmax(freqLength)];
    
               if (freq_rangeNew{j}(2)-freq_rangeNew{j}(1))<=4*freqRes
                   freq_rangeNew{j} = freq_range{j};
               end
                    
               fooof_ind = fooof(freqVals,SpecPower(ind,:)',freq_rangeNew{j},settings,true);
                if knee_flag
                    [exponent(ind,j), offset(ind,j),knee(ind,j),~,~,~,~,r_SQ(ind,j)] = fooof_get_params(fooof_ind,SpecPower(ind,:)',settings);
                else
                    [exponent(ind,j), offset(ind,j),~,~,~,~,~,r_SQ(ind,j)] = fooof_get_params(fooof_ind,SpecPower(ind,:)',settings);
                end
                
                fooof_resultSingElec = [fooof_resultSingElec,fooof_ind];
             end
    %         log_freqs='True';
    %         fooof_plot(fooof_ind);%,log_freqs);
    %         %title(sprintf('elecnum %i ',elecOccipitalNumNew(ind)))
    %         %saveas(gcf,sprintf('case/BLpower_cases %i.jpg',ind))
    %         %plot(fooof_ind.freqs,fooof_ind.power_spectrum)
    %         fooof_ind.peak_params
        fooof_results{ind} = fooof_resultSingElec;

         else
                fooof_results{ind} = getSlopesPSDBaseline_v2(log10(SpecPower(ind,:)),freqVals,centreFreq,freqRangeWidth*ones(1,length(centreFreq)),freqsToAvoid);
                offset(ind,:) = log10(cell2mat(fooof_results{ind}(1,:)));
                exponent(ind,:) = cell2mat(fooof_results{ind}(2,:));
            end

        icheck{icount} = [ind,icount];
        icount=icount+1;
     
         end
     end
     fooof_resultNumVsElecNum=icheck;
     saveFolder = fullfile(folderSourceString,projectName,'FOOOF'); makeDirectory(saveFolder);
     if fooofFlag
         if knee_flag
             saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType 'withRandomTrials' numTrialsForRand 'withKnee.mat']);
         else
             saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType 'withRandomTrials' numTrialsForRand 'withoutKnee.mat']);
         end

     else
         saveFileName = fullfile(saveFolder,[subjectName '_freqWidth_' num2str(freqRangeWidth) '_' refType '_' powerType 'withRandomTrials' numTrialsForRand 'withoutFooof.mat']);
     end
            save(saveFileName,'freqVals','SpecPower','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum','freqRangeWidth','centreFreq','freq_range','settings');


 %fooof_values = struct('freq_range',freq_range,'fooof_results',fooof_results,'exponent',exponent,'offset',offset,'knee',knee,'fooof_resultNumVsElecNum',icheck);
% fooof_values= struct();
% fooof_values.subjectname = subjectName;
% fooof_values.refType = refType;
% fooof_values.freq_range = freq_range;
% fooof_values.fooof_results = fooof_results;
% fooof_values.exponent = exponent;
% fooof_values.offset= offset;
% fooof_values.knee = knee;
% fooof_values.r_square = r_SQ;
% fooof_values.fooof_resultNumVsElecNum = icheck;

 if fooofPlotFlag
     %close all
     figLocation = fullfile(foldername,'FOOOF','plots'); makeDirectory(figLocation);
     figure(1);clf;
        for ind = 1:64
                subplot(8,8,ind);
                 if ~isempty(fooof_results{ind})
                     for jFreq = 1:length(freq_range)
                        fooof_plot(fooof_results{ind}(jFreq),false);
                     end
                   %xlabel(num2str(r_SQ(icheck{ind}(1))))
                    
                 end
                 title(ind);
        end
        sgtitle([subjectName 'freqRangeWidth: ' num2str(freqRangeWidth)])
    set(gcf, 'Position', get(0, 'Screensize')); 
     saveas(gcf,fullfile(figLocation,[subjectName 'freqRangeWidth_ ' num2str(freqRangeWidth) '_1-64_' refType '_' powerType '.fig']));
     
     if length(fooof_results)>64
         figure(2);clf;
             for ind = 1:48
                subplot(8,6,ind);
                 if ~isempty(fooof_results{ind+64})
                     for jFreq = 1:length(freq_range)
                        fooof_plot(fooof_results{ind+64}(jFreq),false);
                     end
                   %xlabel(num2str(r_SQ(icheck{ind}(1))))
                    
                 end
                 title(ind+64);
             end
             sgtitle([subjectName ' 65-112 freqRangeWidth: ' num2str(freqRangeWidth)])
            set(gcf, 'Position', get(0, 'Screensize')); 
         saveas(gcf,fullfile(figLocation,[subjectName  'freqRangeWidth_' num2str(freqRangeWidth) '_64-112_' refType '_' powerType '.fig']));
     end
 end

    if rmsPlotFlag    r_SQFull{iFreqRangeWidth} = r_SQ;     end
end

if rmsPlotFlag
    figure(3);
    xbins = 0.05:.1:.95;
    clf;for i=1:length(centreFreq)
    subplot(4,4,i)
    r_SQPlot = [];
    for iFreqRange = 1:length(allFreqRangeWidths)
        r_SQPlot = [r_SQPlot , r_SQFull{iFreqRange}(:,i)];
    end
    %hist(r_SQFull{iFreqRange}(:,i),xbins); hold on;
    hist(r_SQPlot,xbins); hold on;
    title(centreFreq(i));
    xlim([0.4 1]);
    end
    %legend(mat2str(allFreqRangeWidths));
    sgtitle( [subjectName ' R^2 freqRangeWidths: ' mat2str(allFreqRangeWidths)]);
     figLocation = fullfile(foldername,'FOOOF','rmsPlots'); makeDirectory(figLocation);
     set(gcf, 'Position', get(0, 'Screensize')); 
       saveas(gcf,fullfile(figLocation,[subjectName '_' refType '_' powerType '.fig']));
end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%rough%%%%%%%%%%%%%%%%%%%%
% figure();
%  for i=1:15
% subplot(4,4,i)
% histogram(r_SQ(:,i),'NumBins',10,'BinLimits',[0,1]);
% end
% end
% 
% exponent2 = exponent;
% r_SQ2 = r_SQ;
% fooof_results2 = fooof_results;
% 
%  figure(1); clf;
%         for ind = 1:64
%                 subplot(8,8,ind);
%                  if ~isempty(fooof_results{ind})
%                      for jFreq = 1:length(freq_range)
%                         fooof_plot(fooof_results{ind}(jFreq),false);
%                      end
%                    %xlabel(num2str(r_SQ(icheck{ind}(1))))
%                     
%                  end
%                  title(ind);
%         end
% sgtitle( [subjectName ' 3 tapers freqWidth:' num2str(freqRangeWidth)]);
% 
% figure(2);
% clf;
% for jFreq = 1:length(freq_range)
%                         fooof_plot(fooof_results{24}(jFreq),false);
%                      end
% title('3 tapers');
% 
% figure(3);
% clf;
% plot(centreFreq,r_SQ(24,:),'lineWidth',1)
% hold on;
% hold on;plot(centreFreq,exponent(24,:),'lineWidth',1)
% plot(centreFreq,0.9*(centreFreq.^0),'--k','lineWidth',1);
% plot(centreFreq,(centreFreq.^0),'--k','lineWidth',1);
% plot(centreFreq,0.8*(centreFreq.^0),'--k','lineWidth',1);
% legend('R^2' ,'exponent');
% xlim([10 110]);
% 
% figure(4);
% xbins = 0.05:.1:.95;
% clf;for i=1:16
% subplot(4,4,i)
% hist([r_SQ1(:,i), r_SQ2(:,i), r_SQ3(:,i)],xbins);
% title(freq_range{i});
% xlim([0.5 1]);
% end
% legend('1 taper','2 tapers', '3 tapers');
% sgtitle( [subjectName ' R^2 freqRangeWidth:' num2str(freqRangeWidth)]);
% folderSave = 'C:\Users\Srishty\OneDrive - Indian Institute of Science\Documents\supratim\TLS\TLSAEEGProjectPrograms-master\ageProjectCodes\fooof_AgeProjectCodes\prob n sol\Apt Rsq and Tapers';
% save(fullfile(folderSave,[subjectName '_freqWidth_' num2str(freqRangeWidth)] ),'exponent1','exponent2','exponent3','r_SQ1','r_SQ2','r_SQ3','fooof_results1','fooof_results2','fooof_results3');
% 
