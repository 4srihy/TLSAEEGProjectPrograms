%To get the fooof plot for individual data already stored in TLSAEEG
%analysed data folder


%%%%%%%%getting subject details%%%%%%%%%%%

function saveFooofParamsIndividualSubjects_v2(folderSourceString,subjectName,projectName,refType,protocolType,freq_range, powerType, knee_flag,fooofPlotFlag,specPowerFlag)

if ~exist('knee_flag','var')    knee_flag=1;            end
if ~exist('freq_range','var')   freq_range = [5,70];    end
if ~exist('fooofPlotFlag','var')    fooofPlotFlag=0;    end
if ~exist('powerType','var')    powerType = 'BL';       end
if ~exist('specPowerFlag','var')    specPowerFlag=0;    end

foldername = fullfile(folderSourceString,projectName);
load(fullfile(foldername,protocolType,[subjectName '_' refType '_stRange_250_750.mat']));

if strcmp(powerType,'BL')
    SpecPower= nanmean(blPowerVsFreqTopo,3);
elseif strcmp(powerType,'ST')
    SpecPower = nanmean(stPowerVsFreqTopo,3);
else
    disp("Specify the power type");
end
%for g=[43]   SpecPower(g,:) = nan; end
settings = struct();
settings.peak_width_limits = [4,8];
settings.max_n_peaks = 5;
if knee_flag      settings.aperiodic_mode = 'knee';      end
%settings.peak_threshold = 0.5;
fooof_results=[];
exponent = []; offset=[]; knee=[]; r_SQ=[]; icheck=[];
icount=1;
for ind = 1:size(SpecPower,1)   %index of the cases/control
    if isnan(SpecPower(ind,:))
        exponent(ind) = nan; offset(ind)=nan; knee(ind)=nan; r_SQ(ind)=nan;
    else
        fooof_ind = fooof(freqVals,SpecPower(ind,:)',freq_range,settings,true);
        if knee_flag
            [exponent(ind), offset(ind),knee(ind),~,~,~,~,r_SQ(ind)] = fooof_get_params(fooof_ind,SpecPower(ind,:)',settings);
        else
            [exponent(ind), offset(ind),~,~,~,~,~,r_SQ(ind)] = fooof_get_params(fooof_ind,SpecPower(ind,:)',settings);
        end
        fooof_results = [fooof_results,fooof_ind];
        log_freqs='True';
        %         fooof_plot(fooof_ind);%,log_freqs);
        %         %title(sprintf('elecnum %i ',elecOccipitalNumNew(ind)))
        %         %saveas(gcf,sprintf('case/BLpower_cases %i.jpg',ind))
        %         %plot(fooof_ind.freqs,fooof_ind.power_spectrum)
        %         fooof_ind.peak_params
        icheck{icount} = [ind,icount];
        icount=icount+1;

    end
end
fooof_resultNumVsElecNum=icheck;

if knee_flag
    if specPowerFlag
        save(fullfile(folderSourceString,projectName,'FOOOF',[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '_' powerType '.mat']),'freq_range','SpecPower','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum');
    else
        save(fullfile(folderSourceString,projectName,'FOOOF',[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '_' powerType '.mat']),'freq_range','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum');
    end
else
    if specPowerFlag
        save(fullfile(folderSourceString,projectName,'FOOOF',[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '_' powerType 'WithoutKnee.mat']),'freq_range','SpecPower','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum');
    else
        save(fullfile(folderSourceString,projectName,'FOOOF',[subjectName '_freqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_' refType '_' powerType 'WithoutKnee.mat']),'freq_range','fooof_results','exponent','offset','knee','r_SQ','fooof_resultNumVsElecNum');
    end
end

if fooofPlotFlag
    close all
    figLocation = fullfile(foldername,'FOOOF','plots');
    figure(1);
    for i=1:8
        for j=1:8
            ind = 8*(i-1)+j;
            subplot(8,8,ind);
            if (ind) < length(fooof_results)
                fooof_plot(fooof_results(ind),true);
                xlabel(num2str(r_SQ(icheck{ind}(1))))
                title(icheck{ind}(1))
            end
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf,fullfile(figLocation,[subjectName '_FreqRange_'  num2str(freq_range(1)) '_' num2str(freq_range(2)) '_1-64_' refType '_' powerType '.png']));

    if length(fooof_results)>64
        figure(2);
        for i=1:8
            for j=1:8
                ind = 8*(i-1)+j;
                subplot(8,8,ind);
                if (ind+64) < length(fooof_results)
                    fooof_plot(fooof_results(ind+64),true);
                end
            end
        end
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,fullfile(figLocation,[subjectName '_FreqRange_' num2str(freq_range(1)) '_' num2str(freq_range(2)) '_64-112_' refType '_' powerType '.png']));
    end
end
end

