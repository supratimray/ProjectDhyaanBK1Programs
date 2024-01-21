% plotFigures
comparisonStr = 'paired';
protocolNames = {'G2'}; refChoices = {'G2'};
% protocolNames = [{'G1'} {'G2'}]; refChoices = [{'G1'} {'G2'}];

analysisChoice = 'st';

badEyeCondition = 'ep';
badTrialVersion = 'v8';
badElectrodeRejectionFlag = 1;

stRange = [0.25 1.25]; % hard coded for now

freqRangeList{1} = [8 13];
freqRangeList{2} = [24 34];
freqRangeList{3} = [35 65];

axisRangeList{1} = [0 100]; axisRangeName{1} = 'Freq Lims (Hz)';
axisRangeList{2} = [-2.5 2.5]; axisRangeName{2} = 'YLims';
axisRangeList{3} = [-1.5 1.5]; axisRangeName{3} = 'cLims (topo)';

cutoffList = [3 30];
useMedianFlag = 0;

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
hAllPlots = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
groupPos = 1; % Occipital
freqPos = 2; % Slow gamma

numProtocols = length(protocolNames);
if numProtocols==1
    [psdDataToReturn,powerDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] = displayPowerDataAllSubjects(subjectNameLists,protocolNames{1},analysisChoice,refChoices{1},badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,0);
    logPSDData = psdDataToReturn{groupPos};
    logPower = powerDataToReturn{groupPos,freqPos};

else % Combine
    logPSDDataTMP = cell(1,numProtocols);
    logPowerTMP = cell(1,numProtocols);
    goodSubjectNameListsTMP = cell(1,numProtocols);
    for i=1:numProtocols
        [psdDataTMP,powerDataTMP,goodSubjectNameListsTMP{i},topoplotDataTMP,freqVals] = displayPowerDataAllSubjects(subjectNameLists,protocolNames{i},analysisChoice,refChoices{i},badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,0);
        logPSDDataTMP{i} = psdDataTMP{groupPos};
        logPowerTMP{i} = powerDataTMP{groupPos,freqPos};
    end

    % Average data across protocols
    logPSDData = cell(1,2);
    logPower = cell(1,2);

    for i=1:2 % meditator/control

        % Get subjects for each protocol
        subjectNamesTMP = cell(1,numProtocols);
        for j=1:numProtocols
            subjectNamesTMP{j} = goodSubjectNameListsTMP{j}{groupPos,i};
        end
        
        % Get common subjects
        commonSubjects = subjectNamesTMP{1};
        for j=2:numProtocols
            commonSubjects = intersect(commonSubjects,subjectNamesTMP{j},'stable');
        end

        % Generate average data across protocols for each common subject
        numCommonSubjects = length(commonSubjects);

        logPSDCommon = zeros(numCommonSubjects,length(freqVals));
        logPowerCommon = zeros(1,numCommonSubjects);
        for j=1:numCommonSubjects
            name = commonSubjects{j};

            psdTMP=[]; powerTMP=[];
            for k=1:numProtocols
                pos = find(strcmp(name,subjectNamesTMP{k}));
                psdTMP = cat(3,psdTMP,logPSDDataTMP{k}{i}(pos,:));
                powerTMP = cat(2,powerTMP,logPowerTMP{k}{i}(pos));
            end

            logPSDCommon(j,:) = squeeze(mean(psdTMP,3));
            logPowerCommon(j) = mean(powerTMP);
        end

        logPSDData{i} = logPSDCommon;
        logPower{i} = logPowerCommon;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [1 0 0];
displaySettings.colorNames(2,:) = [0 1 0];
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

freqLims = axisRangeList{1};
yLimsPSD = axisRangeList{2};
cLimsTopo = axisRangeList{3};

hPSD = subplot(221);
displayAndcompareData(hPSD,logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,~pairedDataFlag);
xlim(hPSD,freqLims);

% Add lines in PSD plots
for k=1:2
    line([freqRangeList{freqPos}(k) freqRangeList{freqPos}(k)],yLimsPSD,'color','k','parent',hPSD);
end

if ~strcmp(refChoices{1},'none')
    line([0 freqVals(end)],[0 0],'color','k','parent',hPSD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPower = subplot(223);
displaySettings.plotAxes = hPower;
if ~useMedianFlag
    displaySettings.parametricTest = 1;
else
    displaySettings.parametricTest = 0;
end

displaySettings.showYTicks=1;
displaySettings.showXTicks=1;
displayViolinPlot(logPower,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);