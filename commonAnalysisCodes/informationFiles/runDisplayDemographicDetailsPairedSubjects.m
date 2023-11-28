%%Code for generating paired dots with line
projectName = 'BK1';
[subjectNameList,expDateList,labelList,ageList,genderList,educationList,mcList] = getDemographicDetails(projectName);

pairedSubjectNameList = getPairedSubjectsBK1;
numPairs = size(pairedSubjectNameList,1);
numGroup = size(pairedSubjectNameList,2);
strList = [{'M'} {'C'}];
d = struct();
param = [{'age'} {'education'} {'mc'} {'expDate'}];
paramUnit = [{' (Years)'} {' (Years)'} {' (Days)'} {''}];
numParam = length(param);

for iGroup = 1:numGroup
    for iPair=1:numPairs
        pos = strcmp(pairedSubjectNameList{iPair,iGroup},subjectNameList);
        d.age{iGroup}(iPair) = ageList(pos);
        d.education{iGroup}(iPair) = educationList(pos);
        d.mc{iGroup}(iPair) = mcList(pos);
        x = expDateList{pos};
        d.expDate{iGroup}(iPair)= datenum([x(1:2) '/' x(3:4) '/' x(5:6)],'dd/mm/yy');%datetime(x,'inputFormat','ddmmyy');%
    end
end

figure;
numRows = 2;    numColumns = numParam;
for iParam = 1:length(param)
    %% data Plot
    data = d.(param{iParam});
    dataDiff = data{1}-data{2};
    if iParam ==4
        for i=1:numGroup
            data{i} = datetime(data{i},'ConvertFrom','datenum');
        end
    end

    hPlot(1) = subplot(numRows,numColumns,iParam);
    plotMatchedPairedDots(hPlot(1),data);
    xlim([0 3]);
    set(hPlot(1),'XTick',[1 2],'XTickLabel',[{strList{1}} {strList{2}}],'TickLength',[0.04 0.02],'LineWidth',1,'FontWeight','bold');
    ylabel([upper(param{iParam}(1)) param{iParam}(2:end) paramUnit{iParam}],FontSize=12,FontWeight="bold");
    title([upper(param{iParam}(1)) param{iParam}(2:end)],FontSize=14,FontWeight="bold");

    %% diff plot
    hPlot(2) = subplot(numRows,numColumns,iParam+numColumns);

    [~,p] = ttest(dataDiff);
    err = nanstd(dataDiff)./sqrt(numPairs);

    swarmchart(hPlot(2), ones(1,numPairs),dataDiff,20,[1 0.5 0],'filled','MarkerFaceAlpha',0.7);
    hold on;
    xlim([0 2]);

    plot(hPlot(2).XLim,[0 0],'k','LineWidth',1);
    plot(1,nanmean(dataDiff),'.k',MarkerSize=18);
    er = errorbar(1,nanmean(dataDiff),err,err,'LineWidth',1.5,'Color','k');
    ylim1  = max(abs(hPlot(2).YLim));
    ylim([-ylim1 ylim1]);

    % boxchart(ones(1,numPairs),dataDiff)
    pPosition=ylim1;
    if p<0.001
        text(0.5,pPosition,['p = ' num2str(p,'%.1e')]);
    else
        text(0.5,pPosition,['p = ' num2str(p,'%.3f')]);
    end

    set(hPlot(2),'XTickLabel',[],'FontWeight','bold','LineWidth',1.5);
    if iParam~=4
        ylabel([upper(param{iParam}(1)) param{iParam}(2:end) ' Difference' paramUnit{iParam}],FontSize=12,FontWeight="bold");
    else
        ylabel([upper(param{iParam}(1)) param{iParam}(2:end) ' Difference' paramUnit{iParam-1}],FontSize=12,FontWeight="bold");
    end

end

function plotMatchedPairedDots(hPlot,data)
numGroups = length(data);
numPairs = length(data{1});
for i=1:numGroups
    plot(hPlot,i*ones(1,numPairs),data{i},'o','MarkerFaceColor','b','MarkerEdgeColor','k');
    hold on;
end

for u= 1:length(data{1})
    plot(hPlot,[1 2],[data{1}(u) data{2}(u)],'k','LineWidth',0.8);
    hold on;
end
end