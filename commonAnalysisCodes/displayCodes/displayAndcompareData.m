function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,nonMatchedFlag,compareRefFlag,refIndex,omitFreqFlag,omitFreqRange)

if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end
if ~exist('compareRefFlag','var');          compareRefFlag=0;           end
if ~exist('refIndex','var');                refIndex=1;                 end
if ~exist('omitFreqFlag', 'var');           omitFreqFlag = 0;           end
if ~exist('omitFreqRange', 'var');          omitFreqRange = [47 53];    end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1)));
else
    getLoc = @(g)(squeeze(mean(g,1)));
end

numGroups = length(data);

axes(hPlot);

if omitFreqFlag
    % omit line noise flag
    % If the flag is on, it will not display the data at 'omitFreqRange'
    omitIdx = xs >= omitFreqRange(1) & xs <= omitFreqRange(2);
else
    omitIdx = false(size(xs));
end

for i = 1:numGroups
    clear bootStat mData sData
    mData = getLoc(data{i});
    if useMedianFlag
        rng('default'); % to get the same SEM each time
        bootStat = bootstrp(1000, getLoc, data{i});
        sData = std(bootStat, 0, 1, "omitnan");
    else
        sData = std(data{i}, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data{i}(:,1))));
    end

    % updating the indices as per omitFreqFlag
    xs_valid = xs;
    xs_valid(omitIdx) = NaN;
    mData(omitIdx) = NaN;
    sData(omitIdx) = NaN;

    % plot the data and display significance
    plotData(xs_valid, mData, sData, displaySettings.colorNames(i,:));
end

set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

% Significance testing
if displaySignificanceFlag
    if numGroups > 2
        % Height for significance patches
        yRange = diff(yLims);
        yBase = yLims(1);
        yOffset = yRange/20;

        for i=1:length(xs)
            % Extract data at current frequency point
            freqData = cellfun(@(x) x(:,i), data, 'UniformOutput', false);

            % Get x coordinates for patch
            xMidPos = xs(i);
            if i==1
                xBegPos = xMidPos;
            else
                xBegPos = xMidPos-(xs(i)-xs(i-1))/2;
            end
            if i==length(xs)
                xEndPos = xMidPos;
            else
                xEndPos = xMidPos+(xs(i+1)-xs(i))/2;
            end
            xVals = [xBegPos xEndPos xEndPos xBegPos]';

            if compareRefFlag
                % Perform t-tests against reference group for other groups
                patchCount = 0;
                for g = 1:numGroups
                    if g ~= refIndex
                        if useMedianFlag
                            p = ranksum(freqData{refIndex}, freqData{g});
                        else
                            if nonMatchedFlag
                                [~,p] = ttest2(freqData{refIndex}, freqData{g});
                            else
                                [~,p] = ttest(freqData{refIndex}, freqData{g});
                            end
                        end

                        % Stack patches consecutively
                        yVals = yBase + [0 0 yOffset yOffset] + patchCount*yOffset;
                        patchCount = patchCount + 1;

                        % Plot significance patches
                        if p < 0.01
                            patch(xVals, yVals, 'k', 'linestyle', 'none');
                        elseif p < 0.05
                            patch(xVals, yVals, 'c', 'linestyle', 'none');
                        end
                    end
                end
            else
                % Perform ANOVA/Kruskal-Wallis across all groups
                allData = [];
                allIDs = [];
                for j=1:numGroups
                    allData = cat(1, allData, freqData{j});
                    allIDs = cat(1, allIDs, j*ones(size(freqData{j})));
                end

                if useMedianFlag
                    p = kruskalwallis(allData, allIDs, 'off');
                else
                    [p,~] = anova1(allData, allIDs, 'off');
                end

                % Single row of significance patches for overall comparison
                yVals = yBase + [0 0 yOffset yOffset];
                if p < 0.01
                    patch(xVals, yVals, 'k', 'linestyle', 'none');
                elseif p < 0.05
                    patch(xVals, yVals, 'c', 'linestyle', 'none');
                end
            end
        end
    else
        allData = [];
        allIDs = [];
        for j = 1:numGroups
            allData = cat(1, allData, data{j});
            allIDs = cat(1, allIDs, j + zeros(size(data{j},1), 1));
        end

        for i = 1:length(xs)
            if omitIdx(i), continue; end  % Skipping omitted frequencies

            if useMedianFlag
                if nonMatchedFlag % unpaired
                    p = kruskalwallis(allData(:,i), allIDs, 'off');
                else  % paired
                    group1 = allData(allIDs == 1, i);
                    group2 = allData(allIDs == 2, i);
                    p = signrank(group1, group2);
                end
            else  % mean testing
                if nonMatchedFlag % unpaired
                    [~, p] = ttest2(data{1}(:,i), data{2}(:,i));
                else  % paired
                    [~, p] = ttest(data{1}(:,i), data{2}(:,i));
                end
            end

            % Getting patch coordinates
            for j = 1:numGroups
                if min(data{j}) < diff(yLims) / 10
                    yVals = yLims(1) + [0 0 diff(yLims)/60 diff(yLims)/60]; % to decrease the height of significance pathes when data is too low
                else
                    yVals = yLims(1) + [0 0 diff(yLims)/20 diff(yLims)/20];
                end
            end
            clear xMidPos xBegPos xEndPos
            xMidPos = xs(i);
            if i == 1
                xBegPos = xMidPos;
            else
                xBegPos = xMidPos - (xs(i) - xs(i-1)) / 2;
            end
            if i == length(xs)
                xEndPos = xMidPos;
            else
                xEndPos = xMidPos + (xs(i+1) - xs(i)) / 2;
            end
            clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';

            % Plotting significance patches only if NOT in omitted range
            if ~omitIdx(i)
                if p < 0.05, patch(xVals, yVals, 'c', 'linestyle', 'none'); end
                if p < 0.01, patch(xVals, yVals, 'k', 'linestyle', 'none'); end
            end
        end
    end
end
end

% Helper function to plot the data and significance. 
function plotData(xs, mData, sData, color)
nanIdx = isnan(xs);
hold on;

% Identifying segments (continuous regions without NaNs)
segmentStart = find(diff([false; ~nanIdx(:); false]) == 1);
segmentEnd = find(diff([false; ~nanIdx(:); false]) == -1) - 1;

for seg = 1:length(segmentStart)
    idxRange = segmentStart(seg):segmentEnd(seg);

    xVals = xs(idxRange)';
    yLower = (mData(idxRange) - sData(idxRange))';
    yUpper = (mData(idxRange) + sData(idxRange))';

    % Plotting main curve
    plot(xVals, mData(idxRange), 'color', color, 'linewidth', 1);

    % Using column vectors for patch function
    patch([xVals; flipud(xVals)], [yLower; flipud(yUpper)], ...
        color, 'linestyle', 'none', 'FaceAlpha', 0.4);
end
end