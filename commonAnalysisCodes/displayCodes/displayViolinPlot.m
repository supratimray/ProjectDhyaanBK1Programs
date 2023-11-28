function figHandles=displayViolinPlot(dataArray,colorArray,showData,plotQuartiles,showSignificance,pairedDataFlag)
% uses statistical toolbox function ksdensity
% adapted from "https://github.com/bastibe/Violinplot-Matlab"
% Ref for violinplot: "stat.cmu.edu/~rnugent/PCMI2016/papers/ViolinPlots.pdf"

% get the data in the cell array
Y{:,1}=dataArray{1,1};
Y{:,2}=dataArray{1,2};
bandwidth = [];
alpha = 0.3;
dataMarkerSize = 24;
MedianMarkerSize = 30;

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
    figHandles{pos} = patch([pos+density*width pos-density(end:-1:1)*width], ...
        [value value(end:-1:1)],...
        colorArray{pos},'FaceAlpha',alpha);
    hold on

    if showData
        [~, unique_idx] = unique(value);
        jitterstrength = interp1(value(unique_idx), density(unique_idx)*width, data, 'linear','extrap');
        jitter = 2*(rand(size(data))-0.5);
        xPosData = pos + jitter.*jitterstrength;
        if pairedDataFlag xPosDataGroups(pos,:) = xPosData; end
        scatter(xPosData, data, dataMarkerSize, 'filled','MarkerFaceColor',colorArray{pos});
    end

    if plotQuartiles
        quartiles = quantile(data, [0.25, 0.5, 0.75]);
        IQR = quartiles(3) - quartiles(1);
        lowhisker = quartiles(1) - 1.5*IQR;
        lowhisker = max(lowhisker, min(data(data > lowhisker)));
        hiwhisker = quartiles(3) + 1.5*IQR;
        hiwhisker = min(hiwhisker, max(data(data < hiwhisker)));

        patch(pos+[-1,1,1,-1]*BoxWidth, ...
            [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], ...
            [0.5 0.5 0.5]);
        plot([pos pos], [lowhisker hiwhisker]);
        scatter(pos, quartiles(2), MedianMarkerSize, [1 1 1], 'filled');
    end
end

if pairedDataFlag
    for i=1:length(xPosDataGroups)
        xPosLine = xPosDataGroups(:,i)';
        yPosLine = [Y{1,1}(i) Y{1,2}(i)];
        plot(xPosLine,yPosLine,'Color',[0.5 0.5 0.5]);
    end
end

if showSignificance
    if pairedDataFlag
        % [h, p, ci, stats] = ttest(Y{:,1},Y{:,2});
        % Perform Wilcoxon signed-rank test
        [p, h, stats] = signrank(Y{:,1},Y{:,2});
    else
        % [h, p, ci, stats] = ttest2(Y{:,1},Y{:,2});
        % Perform Mann-Whitney U test
        [p,h,stats] = ranksum(Y{:,1},Y{:,2});
    end
    xPos = [1:size(Y,2)];
    commonMax=max([max(Y{:,1}) max(Y{:,2})]);
    commonMin = min([max(Y{:,1}) min(Y{:,2})]);
    yPos = [commonMax+1 commonMax+1];
    plot(xPos,yPos);
    set(gca,'YLim',[commonMin-3 commonMax+4]);

    if p<0.05
        text(mean(xPos),commonMax+1.38,'*','FontSize',20);
        text(mean(xPos)+0.01,commonMax+1.8,[' (' num2str(round(p,4)) ')'],'FontSize',16);
    else
        text(mean(xPos),commonMax+1.65,'N.S.','FontSize',16);
        text(mean(xPos)+0.06,commonMax+1.65,[' (' num2str(round(p,4)) ')'],'FontSize',16);
    end
end
end