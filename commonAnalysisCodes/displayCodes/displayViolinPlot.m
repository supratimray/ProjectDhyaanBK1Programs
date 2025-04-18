function ax=displayViolinPlot(dataArray,colorArray,showData,plotCentralTendency,showSignificance,pairedDataFlag,displaySettings)
% displayViolinPlot makes the violin plots using the kernel desnsity estimate
% of the original data, kernel density is estimated using the Matlab statistical toolbox function, ksdensity
% The code is adapted from the original source, "https://github.com/bastibe/Violinplot-Matlab"
% To know about the violinplot please refer to "stat.cmu.edu/~rnugent/PCMI2016/papers/ViolinPlots.pdf"

% Parametric testing is now enabled by default.
% Improved spacing for displaying p-values in ANOVA.

if ~exist('plotCentralTendency','var');           plotCentralTendency=0;                          end
if ~exist('displaySettings','var');               displaySettings=struct();                       end
if ~isfield(displaySettings,'alpha');             displaySettings.alpha=0.3;                      end
if ~isfield(displaySettings,'dataMarkerSize');    displaySettings.dataMarkerSize=12;              end
if ~isfield(displaySettings,'medianMarkerSize');  displaySettings.medianMarkerSize=20;            end
if ~isfield(displaySettings,'textFontSize');      displaySettings.textFontSize=8;                 end
if ~isfield(displaySettings,'yPositionLine');     displaySettings.yPositionLine=0.5;              end
if ~isfield(displaySettings,'xPositionText');     displaySettings.xPositionText=0.5;              end
if ~isfield(displaySettings,'plotAxes');          displaySettings.plotAxes=gca;                   end
if ~isfield(displaySettings,'showYTicks');        displaySettings.showYTicks=0;                   end
if ~isfield(displaySettings,'setYLim');           displaySettings.setYLim=[-7 7];                 end
if ~isfield(displaySettings,'commonYLim');        displaySettings.commonYLim=0;                   end
if ~isfield(displaySettings,'showXTicks');        displaySettings.showXTicks=1;                   end
if ~isfield(displaySettings,'xTickLabels');       displaySettings.xTickLabels=[{'Med'},{'Con'},{'Ref'}];  end
if ~isfield(displaySettings,'parametricTest');    displaySettings.parametricTest=1;               end
if ~isfield(displaySettings,'tickLengthMedium');  displaySettings.tickLengthMedium=[0.025 0];     end
if ~isfield(displaySettings,'plotQuartiles');     displaySettings.plotQuartiles=0;                end
if ~isfield(displaySettings,'BoxWidth');          displaySettings.BoxWidth=0.005;                 end
if ~isfield(displaySettings,'scaleFactor');       displaySettings.scaleFactor=2;                  end
if ~isfield(displaySettings,'medianFlag');        displaySettings.medianFlag=0;                   end

alpha            = displaySettings.alpha;
dataMarkerSize   = displaySettings.dataMarkerSize ;
medianMarkerSize = displaySettings.medianMarkerSize;
textFontSize     = displaySettings.textFontSize;
yPositionLine    = displaySettings.yPositionLine;
xPositionText    = displaySettings.xPositionText;
setYLim          = displaySettings.setYLim;
commonYLim       = displaySettings.commonYLim;
showXTicks       = displaySettings.showXTicks;
xTickLabels      = displaySettings.xTickLabels;
parametricTest   = displaySettings.parametricTest;
plotQuartiles    = displaySettings.plotQuartiles;
BoxWidth         = displaySettings.BoxWidth;
scaleFactor      = displaySettings.scaleFactor;
useMedianFlag    = displaySettings.medianFlag;

% get the data in the cell array
bandwidth = [];
numGroups = size(dataArray, 2);
Y = cell(1,numGroups);
for i = 1:numGroups
    Y{i} = dataArray{1,i};
end

% set Plot Options:
axes(displaySettings.plotAxes);
ax=gca;
set(ax,'TickDir','out','TickLength',displaySettings.tickLengthMedium);
ax.XTick=1:numGroups;
if showXTicks
    ax.XTickLabel = xTickLabels;
else
    ax.XTickLabel = [];
end

% calculate kernel density
centralTendency = zeros(1,numGroups);
semData         = zeros(1,numGroups);
xPosDataGroups  = zeros(numGroups,length(Y{1}));

for pos=1:numGroups
    width = 0.3;
    data = Y{pos};

    % Check if data is all zeros or empty
    if isempty(data) || all(data == 0)
        % Draw a simple line at zero for zero vectors
        plot([pos-width pos+width], [0 0], 'Color', colorArray{pos}, 'LineWidth', 2);
        if showData
            scatter(pos, 0, dataMarkerSize, 'filled', 'MarkerFaceColor', colorArray{pos});
            if pairedDataFlag
                xPosDataGroups(pos,:) = pos; % Set position to exact pos for zero data
            end
        end
        centralTendency(pos) = 0;
        semData(pos) = 0;
        hold on
        continue;
    end

    % Original density calculation for non-zero data
    [density, value] = ksdensity(data,'bandwidth',bandwidth);
    density = density(value >= min(data) & value <= max(data));
    value = value(value >= min(data) & value <= max(data));
    value(1) = min(data);
    value(end) = max(data);
    value = [value(1)*(1-1E-10), value, value(end)*(1+1E-10)];
    density = [0, density, 0];

    % violinWidth and the boxWidth
    width = width/max(density);

    % plot violin plot
    patch([pos+density*width pos-density(end:-1:1)*width], ...
        [value value(end:-1:1)],...
        colorArray{pos},'FaceAlpha',alpha);
    hold on

    if showData
        [~, unique_idx] = unique(value);
        jitterstrength = interp1(value(unique_idx), density(unique_idx)*width, data, 'linear','extrap');
        jitter = 2*(rand(size(data))-0.5);
        xPosData = pos + jitter.*jitterstrength;
        if pairedDataFlag
            xPosDataGroups(pos,:) = xPosData;
        end
        scatter(xPosData, data, dataMarkerSize, 'filled','MarkerFaceColor',colorArray{pos});
    end

    if useMedianFlag
        centralTendency (pos) = median(data,'omitnan');
    else
        centralTendency(pos) = mean(data,'omitnan');
    end
    semData(pos)  = std(data,'omitnan')./sqrt(nnz(~isnan(data)));

    if plotQuartiles
        quartiles = quantile(data, [0.25, 0.5, 0.75]);
        IQR = quartiles(3) - quartiles(1);
        lowhisker = quartiles(1) - 1.5*IQR;
        lowhisker = max(lowhisker, min(data(data > lowhisker)));
        hiwhisker = quartiles(3) + 1.5*IQR;
        hiwhisker = min(hiwhisker, max(data(data < hiwhisker)));

        patch(pos+[-1,1,1,-1]*(BoxWidth+0.005), ...
            [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], ...
            [0 0 0]);
        plot([pos pos], [lowhisker hiwhisker]);
        scatter(pos, quartiles(2), medianMarkerSize, [1 1 1], 'filled');
    end
end

if pairedDataFlag
    for i=1:size(xPosDataGroups,2) % loop through subjects
        for j=1:numGroups-1 % loop through adjacent pairs
            xPosLine = [xPosDataGroups(j,i) xPosDataGroups(j+1,i)];
            yPosLine = [Y{j}(i) Y{j+1}(i)];
            plot(xPosLine,yPosLine,'Color',[0.8 0.8 0.8]);
        end
    end
end

if plotCentralTendency
    for pos=1:numGroups
        patch(pos+[-1,1,1,-1]*(BoxWidth+0.005), ...
            [centralTendency(pos)-semData(pos) centralTendency(pos)-semData(pos) centralTendency(pos)+semData(pos) centralTendency(pos)+semData(pos)], ...
            [1 1 1]);
        scatter(pos, centralTendency(pos), medianMarkerSize+10, [0 0 0], 'filled');
    end
end

if showSignificance
    if numGroups > 2
        % Perform ANOVA/KW test first
        if parametricTest
            [p_anova,~,~] = anova1(cell2mat(Y),[],'off');
        else
            [p_anova,~,~] = kruskalwallis(cell2mat(Y),[],'off');
        end

        % If compareRefFlag is on, also do t-tests against reference
        if isfield(displaySettings, 'compareRefFlag') && displaySettings.compareRefFlag && isfield(displaySettings, 'refIndex')
            pairwise_p = zeros(1, numGroups);
            for g = 1:numGroups
                if g ~= displaySettings.refIndex
                    if parametricTest
                        if pairedDataFlag
                            [~,pairwise_p(g)] = ttest(Y{displaySettings.refIndex}, Y{g});
                        else
                            [~,pairwise_p(g)] = ttest2(Y{displaySettings.refIndex}, Y{g});
                        end
                    else
                        pairwise_p(g) = ranksum(Y{displaySettings.refIndex}, Y{g});
                    end
                end
            end
        end
    else
        if pairedDataFlag
            if parametricTest
                [~, p, ~, ~] = ttest(Y{:,1},Y{:,2});
            else
                % Perform Wilcoxon signed-rank test
                [p, ~, ~] = signrank(Y{:,1},Y{:,2});
            end
        else
            if parametricTest
                [~, p, ~, ~] = ttest2(Y{:,1},Y{:,2});
            else
                % Perform Mann-Whitney U test
                [p,~,~] = ranksum(Y{:,1},Y{:,2});
            end
        end
    end

    xPos = 1:numGroups;
    commonMax = max(cellfun(@max, Y));
    commonMin = min(cellfun(@min, Y));

    if commonYLim
        set(ax,'YLim',setYLim);
    else
        set(ax,'YLim',[commonMin-yPositionLine commonMax+yPositionLine*scaleFactor]);
    end

    % Show significance text
    if numGroups > 2
        % Calculate spacing for significance markers with larger gaps
        ySpacing = yPositionLine;
        baseYPos = commonMax + yPositionLine;

        % Show ANOVA result at the top with extra spacing
        text(mean(xPos), baseYPos, ...
            ['ANOVA: ' getSignificanceString(p_anova)], ...
            'FontSize', textFontSize, ...
            'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center');

        % If compareRefFlag is on, show pairwise comparisons
        if isfield(displaySettings, 'compareRefFlag') && displaySettings.compareRefFlag && isfield(displaySettings, 'refIndex')
            % Position significance stars with larger vertical gaps
            for g = 1:numGroups
                if g ~= displaySettings.refIndex && pairwise_p(g) <= 0.05
                    if g < displaySettings.refIndex
                        xPos = g;
                        yPos = baseYPos - 1.5*ySpacing*(displaySettings.refIndex - g);
                    else
                        xPos = g;
                        yPos = baseYPos - 1.5*ySpacing*(g - displaySettings.refIndex);
                    end

                    % Draw connection line to reference with adjusted position
                    plot([xPos displaySettings.refIndex], [yPos-ySpacing/3 yPos-ySpacing/3], ...
                        'Color', [0.5 0.5 0.5], 'LineStyle', ':');

                    % Add significance star with larger font
                    text(xPos, yPos, getSignificanceString(pairwise_p(g)), ...
                        'FontSize', textFontSize * 1.2, ... % Slightly larger font
                        'FontWeight', 'bold', ...
                        'HorizontalAlignment', 'center');
                end
            end
        end
    else
        if p>0.05
            text(mean(xPos)-xPositionText/scaleFactor,commonMax+yPositionLine,['N.S. (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
        elseif p>0.01
            text(mean(xPos)-xPositionText/scaleFactor,commonMax+yPositionLine,['* (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
        elseif p>0.005
            text(mean(xPos)-xPositionText/scaleFactor,commonMax+yPositionLine,['** (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
        else
            text(mean(xPos)-xPositionText/scaleFactor,commonMax+yPositionLine,['*** (' num2str(round(p,3)) ')'],'FontSize',textFontSize,'FontWeight','bold');
        end
    end
end
end

function str = getSignificanceString(p)
if p > 0.05
    str = 'n.s.';
elseif p > 0.01
    str = '*';
elseif p > 0.005
    str = '**';
else
    str = '***';
end
end