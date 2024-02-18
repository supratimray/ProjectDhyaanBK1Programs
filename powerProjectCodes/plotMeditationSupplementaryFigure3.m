% Plots Supplementary Figure 3 for the BK1 Project
% Shows Correlations between stimulus-induced and meditation-induced gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
clear
comparisonStr = 'paired';
fontSize = 16;
figure1 = figure('WindowState','maximized','Color',[1 1 1]);
if strcmp(comparisonStr,'paired')
    pairedSubjectNameList = getPairedSubjectsBK1;
    subjectNameLists{1} = pairedSubjectNameList(:,1);
    subjectNameLists{2} = pairedSubjectNameList(:,2);
    pairedDataFlag      = 1;
else
    [~, meditatorList, controlList] = getGoodSubjectsBK1;
    subjectNameLists{1} = meditatorList;
    subjectNameLists{2} = controlList;
    pairedDataFlag      = 0;
end

badEyeCondition = 'ep';
badTrialVersion = 'v8';
badElectrodeRejectionFlag = 1;

stRange = [0.25 1.25]; % hard coded for now

axisRangeList{1} = [0 100];
axisRangeList{2} = [-2.5 2.5];
axisRangeList{3} = [-1.5 1.5];

cutoffList = [3 30]; useMedianFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Meditation-induced gamma (M1-EO1)    vs stimulus-induced gamma(G1_ST-G1_BL)
% 2. Meditation-induced gamma (M2-G2(bl)) vs stimulus-induced gamma(M2(st)-M2(bl))

protocolNames  = [{'M1'}       {'M1'}         {'M2'}    {'G1'}     {'M2'}   ];
refChoices     = [{'none'}     {'EO1'}        {'G2'}    {'G1'}     {'M2'}   ] ;
analysisChoice = [{'combined'} {'combined'}   {'bl'}    {'st'}     {'st'}   ];
freqRangeList  = {[30 80],      [30 80]       [30 80]   [24 34]    [24 34]  };
groupPosList   = [   1             1             1         1          1     ];

compareIndexesList   = {[2 4] [3 5]};

numProtocol      = length(protocolNames);
logPowerData     = cell(1,numProtocol);
goodSubjectNames = cell(1,numProtocol);
for i=1:numProtocol
    groupPos               = groupPosList(i);
    freqRangeListTmp{1}    = freqRangeList{i};

    [~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolNames{i},analysisChoice{i},refChoices{i},badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeListTmp,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
    logPowerData{1,i}     = powerDataToReturn{groupPos,1};
    goodSubjectNames{1,i} = goodSubjectNameListsToReturn(groupPos,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get figure handles
hf0=getPlotHandles(1,2,[0.1 0.15 0.8 0.8],0.08,0.05);
groupNames = {'Meditators','Controls'};
colorNames(1,:)  = [0.8 0 0.8];
colorNames(2,:)  = [0.25 0.41 0.88];
yLimArray        = {[-2 8], [-2 4]};
xLimArray        = {[-6 6], [-4 8]};

%%%%%%%%%%%%%%%%%%%%%%%%%%% Stimulus Induced vs Meditation-Induced %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(compareIndexesList)
    compareIndexes = compareIndexesList{i};
    r=zeros(1,2);  p=zeros(1,2);  N=zeros(1,2);
    xData=[];      yData=[];
    yshift=0.4;

    for pPos=1:2 % 1 for meditators, 2 for controls
        data1 = logPowerData{compareIndexes(1)}{pPos};
        subjects1 = goodSubjectNames{compareIndexes(1)}{pPos};
        data2 = logPowerData{compareIndexes(2)}{pPos};
        subjects2 = goodSubjectNames{compareIndexes(2)}{pPos};
        [r(pPos),p(pPos),N(pPos)]=getCorrelation(data1,subjects1,data2,subjects2,hf0(1,i),colorNames(pPos,:));
        xData=cat(2,xData,data1);
        yData=cat(2,yData,data2);
    end

    % set the text position
    ylim(hf0(1,i),yLimArray{i});
    xlim(hf0(1,i),xLimArray{i});
    set(hf0(1,i),"FontSize",fontSize-2);
    for pPos=1:2
        xPos=min(xData(:))+5; yPos=yLimArray{i}(2)-yshift;
        text(xPos,yPos,['r = ' num2str(round(r(pPos),2)) ', p = ' num2str(round(p(pPos),2)) ', n = ' num2str(N(pPos))],'Color',colorNames(pPos,:),'parent',hf0(1,i),'FontSize',fontSize-2,'FontWeight','bold');
        if i==1
            yshift=0.8;
        else
            yshift=0.65;
        end
    end
end

% axis Labels
xlabel(hf0(1,1),{'Meditation-induced gamma power (dB)';'(M1-EO1)'},'FontSize',fontSize,'FontWeight','bold');
ylabel(hf0(1,1),{'Stimulus-induced gamma power (dB)';'(G1(St-Bl))'},'FontSize',fontSize,'FontWeight','bold');

xlabel(hf0(1,2),{'Meditation-induced gamma power (dB)'; '(M2(Bl)-G2(Bl))'},'FontSize',fontSize,'FontWeight','bold');
ylabel(hf0(1,2),{'Stimulus-induced gamma power (dB)';'(M2(St-Bl))'},'FontSize',fontSize,'FontWeight','bold');

annotation(figure1,'textbox',[0.1437  0.8940  0.0738  0.0381],'Color',colorNames(1,:),...
    'String',{'Meditators'},...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',fontSize);

annotation(figure1,'textbox',[0.1438  0.8638  0.0737  0.0381],'Color',colorNames(2,:),...
    'String','Controls',...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',fontSize);

% figure captions
annotation(figure1,'textbox',[0.0581  0.9520  0.0738  0.0381],'String','A','LineStyle','none','FontWeight','bold','FontSize',fontSize+4,...
    'FitBoxToText','off');

annotation(figure1,'textbox',[0.4981  0.9520  0.0738  0.0381],'String','B','LineStyle','none','FontWeight','bold','FontSize',fontSize+4,...
    'FitBoxToText','off');

% common change across figure!
set(findobj(figure1,'type','axes'),'box','off'...
    ,'FontWeight','Bold'...
    ,'TickDir','out'...
    ,'TickLength',[0.01 0.01]...
    ,'linewidth',1.2...
    ,'xcolor',[0 0 0]...
    ,'ycolor',[0 0 0]...
    );