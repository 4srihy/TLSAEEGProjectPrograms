function ax=displayViolinPlot_Sri(dataArray,colorArray,showData,plotQuartiles,showSignificance,pairedDataFlag,displaySettings)
% uses statistical toolbox function ksdensity
% adapted from "https://github.com/bastibe/Violinplot-Matlab"
% Ref for violinplot: "stat.cmu.edu/~rnugent/PCMI2016/papers/ViolinPlots.pdf"

if ~exist('displaySettings','var');               displaySettings=struct();                     end
if ~isfield(displaySettings,'alpha');             displaySettings.alpha=0.3;                    end
if ~isfield(displaySettings,'dataMarkerSize');    displaySettings.dataMarkerSize=12;            end
if ~isfield(displaySettings,'medianMarkerSize');  displaySettings.medianMarkerSize=20;          end
if ~isfield(displaySettings,'textFontSize');      displaySettings.textFontSize=8;               end
if ~isfield(displaySettings,'yPositionLine');     displaySettings.yPositionLine=0.5;            end
if ~isfield(displaySettings,'xPositionText');     displaySettings.xPositionText=0.5;              end
if ~isfield(displaySettings,'plotAxes');          displaySettings.plotAxes=gca;                 end
if ~isfield(displaySettings,'showYTicks');        displaySettings.showYTicks=0;                 end
if ~isfield(displaySettings,'setYLim');           displaySettings.setYLim=[-7 7];               end
if ~isfield(displaySettings,'commonYLim');        displaySettings.commonYLim=0;                 end
if ~isfield(displaySettings,'showXTicks');        displaySettings.showXTicks=1;                 end
if ~isfield(displaySettings,'xTickLabels');       displaySettings.xTickLabels=[{'Med'},{'Con'}];    end
if ~isfield(displaySettings,'parametricTest');    displaySettings.parametricTest=0;             end
if ~isfield(displaySettings,'tickLengthMedium');  displaySettings.tickLengthMedium=[0.025 0];   end

% get the data in the cell array
Y{:,1}=dataArray{1,1};
Y{:,2}=dataArray{1,2};
bandwidth = [];
alpha            = displaySettings.alpha;
dataMarkerSize   = displaySettings.dataMarkerSize ;
medianMarkerSize = displaySettings.medianMarkerSize;
textFontSize     = displaySettings.textFontSize;
yPositionLine    = displaySettings.yPositionLine;
xPositionText    = displaySettings.xPositionText;
showYTicks       = displaySettings.showYTicks;
setYLim          = displaySettings.setYLim;
commonYLim       = displaySettings.commonYLim;
showXTicks       = displaySettings.showXTicks;
xTickLabels      = displaySettings.xTickLabels;
parametricTest   = displaySettings.parametricTest;

% set Plot Options:
axes(displaySettings.plotAxes);
ax=gca;
set(ax,'TickDir','out','TickLength',displaySettings.tickLengthMedium);
% if ~showYTicks && commonYLim
%     set(ax,'YTickLabel',[]);
% end

ax.XTick=[1 2];
if showXTicks
    ax.XTickLabel = xTickLabels;
else
    ax.XTickLabel = [];
end

% calculate kernel density
xPosDataGroups = zeros(length(Y),length(Y{:,1}));
for pos=1:size(Y,2)
    width = 0.3;
    data = Y{pos};

    [density, value ] = ksdensity(data,'bandwidth',bandwidth);

    density = density(value >= min(data) & value <= max(data));
    value = value(value >= min(data) & value <= max(data));
    value(1) = min(data);
    value(end) = max(data);
    value = [value(1)*(1-1E-10), value, value(end)*(1+1E-10)];
    density = [0, density, 0];

    % violinWidth and the boxWidth
    width = width/max(density);
    BoxWidth = 0.0*width;

    % plot violin plot
    patch([pos+density*width pos-density(end:-1:1)*width], ...
        [value value(end:-1:1)],...
        colorArray{pos},'FaceAlpha',alpha,'EdgeColor',colorArray{pos});
    hold on

    if showData
        [~, unique_idx] = unique(value);
        jitterstrength = interp1(value(unique_idx), density(unique_idx)*width, data, 'linear','extrap');
        jitter = 2*(rand(size(data))-0.5);
        xPosData = pos + jitter.*jitterstrength;
        if pairedDataFlag xPosDataGroups(pos,:) = xPosData; end
        scatter(xPosData, data, dataMarkerSize, 'filled','MarkerFaceColor',colorArray{pos});
    end
end

if pairedDataFlag
    for i=1:length(xPosDataGroups)
        xPosLine = xPosDataGroups(:,i)';
        yPosLine = [Y{1,1}(i) Y{1,2}(i)];
        plot(xPosLine,yPosLine,'Color',[0.75 0.75 0.75]);
    end
end
for pos=1:size(Y,2)
    
    data = Y{pos};

    if parametricTest
         meanData(pos) = nanmean(data);
    else
         meanData(pos) = nanmedian(data);
    end
    semData(pos) = nanstd(data)./sqrt(length(data));

    if plotQuartiles
        quartiles = quantile(data, [0.25, 0.5, 0.75]);
        IQR = quartiles(3) - quartiles(1);
        lowhisker = quartiles(1) - 1.5*IQR;
        lowhisker = max(lowhisker, min(data(data > lowhisker)));
        hiwhisker = quartiles(3) + 1.5*IQR;
        hiwhisker = min(hiwhisker, max(data(data < hiwhisker)));

        line([pos pos], [lowhisker hiwhisker],'Color',[0 0 0]);
        patch(pos+[-1,1,1,-1]*(BoxWidth+0.07), ...
                    [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], ...
                     colorArray{pos});
        patch(pos+[-1,1,1,-1]*(BoxWidth+0.07), ...
            [  meanData(pos) meanData(pos) meanData(pos) meanData(pos)], ...
            [0 0 0]);

%                 patch(pos+[-1,1,1,-1]*(BoxWidth+0.005), ...
%                     [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], ...
%                     [0 0 0]);
%         patch(pos+[-1,1,1,-1]*(BoxWidth+0.01), ...
%             [meanData(pos)-semData(pos) meanData(pos)-semData(pos) meanData(pos)+semData(pos) meanData(pos)+semData(pos) ], ...
%             [0 0 0]);
        
       % meanData(pos) = quartiles(2);
        % scatter(pos, quartiles(2), medianMarkerSize, [0 0 0], 'filled');
    end
end



% for pos=1:size(Y,2)    
%     scatter(pos, meanData(pos), medianMarkerSize, [0 0 0], 'filled'); 
% end
% er = errorbar([1 2], meanData, semData, semData,'LineWidth',1.5);
% er.Color = [0 0 0];
% er.LineStyle = 'None';
% er.MarkerSize = 8

if showSignificance
    if pairedDataFlag
        if parametricTest
            [h, p, ci, stats] = ttest(Y{:,1},Y{:,2});
        else
            % Perform Wilcoxon signed-rank test
            [p, h, stats] = signrank(Y{:,1},Y{:,2});
        end
    else
        if parametricTest
            [h, p, ci, stats] = ttest2(Y{:,1},Y{:,2});
        else
            % Perform Mann-Whitney U test
            [p,h,stats] = ranksum(Y{:,1},Y{:,2});
        end
    end

    xPos = [1:size(Y,2)];
    commonMax=max([max(Y{:,1}) max(Y{:,2})]);
    commonMin = min([max(Y{:,1}) min(Y{:,2})]);
    yPos = [commonMax+yPositionLine commonMax+yPositionLine];
    %plot(xPos,yPos);

    if commonYLim
        set(ax,'YLim',setYLim);
        % set(ax,'YTickLabel',[]);
    else
        set(ax,'YLim',[commonMin-yPositionLine*2 commonMax+yPositionLine*6]);
    end

    if p>0.05
        text(mean(xPos)-xPositionText/4,setYLim(2)+yPositionLine,['N.S. (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
    elseif p>0.01
        text(mean(xPos)-xPositionText/2,setYLim(2)+yPositionLine,['* (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
    elseif p>0.005
        text(mean(xPos)-xPositionText/2,setYLim(2)+yPositionLine,['** (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
    else
        text(mean(xPos)-xPositionText/2,setYLim(2)+yPositionLine,['*** (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
    end

    %     if p<0.05
    %         text(mean(xPos)-xPositionText/4,commonMax+yPositionLine+0.3,'*','FontSize',textFontSize);
    %         text(mean(xPos)-xPositionText/5,commonMax+yPositionLine+0.3,[' (' num2str(round(p,3)) ')'],'FontSize',textFontSize);
    %     else
    %         text(mean(xPos)-xPositionText/2,commonMax+yPositionLine+0.3,'N.S.','FontSize',textFontSize);
    %         text(mean(xPos),commonMax+yPositionLine+0.3,[' (' num2str(round(p,2)) ')'],'FontSize',textFontSize);
    %     end
end


end