function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,nonMatchedFlag)

if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1)));
else
    getLoc = @(g)(squeeze(mean(g,1)));
end

numGroups = length(data);

axes(hPlot);
for i=1:numGroups
    clear bootStat mData sData
    mData = getLoc(data{i});
    if useMedianFlag
        bootStat = bootstrp(1000,getLoc,data{i});
        sData = std(bootStat);
    else
        sData = std(data{i},[],1)/sqrt(size(data{i},1));
    end

    patch([xs';flipud(xs')],[mData'-sData';flipud(mData'+sData')],displaySettings.colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
    hold on;
    plot(xs,mData,'color',displaySettings.colorNames(i,:),'linewidth',1);
end

set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

if displaySignificanceFlag % Do significance Testing

    allData = [];
    allIDs = [];
    for j=1:numGroups
        allData = cat(1,allData,data{j});
        allIDs = cat(1,allIDs,j+zeros(size(data{j},1),1));
    end

    for i=1:length(xs)
        if useMedianFlag
            p=kruskalwallis(allData(:,i),allIDs,'off');
        else
            if nonMatchedFlag
                [~,p]=ttest2(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            else
                [~,p]=ttest(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
            end
        end
        % Get patch coordinates
        yVals = yLims(1)+[0 0 diff(yLims)/20 diff(yLims)/20];

        clear xMidPos xBegPos xEndPos
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
        clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';

        if (p<0.05)
            patch(xVals,yVals,'c','linestyle','none');
        end
        if (p<0.01)
            patch(xVals,yVals,'k','linestyle','none');
        end
    end
end
end