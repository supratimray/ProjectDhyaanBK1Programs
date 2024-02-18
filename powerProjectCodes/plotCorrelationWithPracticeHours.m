%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comparisonStr = 'paired';

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
% Meditation-induced raw broadband gamma
groupPos = 2; % Fronto-central (choose 1 for occipital)
protocolName = {'M1'}; refChoice = {'none'}; analysisChoice = 'combined'; freqRangeList{1} = [30 80]; 
[~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
logPowerMeditationGammaRaw = powerDataToReturn{groupPos,1};
goodSubjectNamesMeditationGammaRaw = goodSubjectNameListsToReturn(groupPos,:);

% Meditation-induced broadband gamma in M1 compared to EO1
groupPos = 1; % Fronto-central (choose 1 for occipital)
protocolName = {'M1'}; refChoice = {'EO1'}; analysisChoice = 'combined'; freqRangeList{1} = [30 80]; 
[~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
logPowerMeditationGamma = powerDataToReturn{groupPos,1};
goodSubjectNamesMeditationGamma = goodSubjectNameListsToReturn(groupPos,:);

% Stimulus-induced narrowband gamma
groupPos = 1; % Occipital
protocolName = {'G1'}; refChoice = {'G1'}; analysisChoice = 'st'; freqRangeList{1} = [24 34]; 
[~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
logPowerStimulusGamma = powerDataToReturn{groupPos,1};
goodSubjectNamesStimulusGamma = goodSubjectNameListsToReturn(groupPos,:);

% Meditation-induced broadband gamma in M2 compared to G2
groupPos = 1; % Fronto-central (choose 1 for occipital)
protocolName = {'M2'}; refChoice = {'G2'}; analysisChoice = 'bl'; freqRangeList{1} = [30 80]; 
[~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
logPowerMeditationGammaM2 = powerDataToReturn{groupPos,1};
goodSubjectNamesMeditationGammaM2 = goodSubjectNameListsToReturn(groupPos,:);

% Stimulus-induced narrowband gamma in M2
groupPos = 1; % Fronto-central (choose 1 for occipital)
protocolName = {'M2'}; refChoice = {'M2'}; analysisChoice = 'st'; freqRangeList{1} = [24 34]; 
[~,powerDataToReturn,goodSubjectNameListsToReturn] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,[],pairedDataFlag,0);
logPowerStimulusGammaM2 = powerDataToReturn{groupPos,1};
goodSubjectNamesStimulusGammaM2 = goodSubjectNameListsToReturn(groupPos,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make violin plots which should match the results in previous plots
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [1 0 0];
displaySettings.colorNames(2,:) = [0 1 0];
if ~useMedianFlag
    displaySettings.parametricTest = 1;
else
    displaySettings.parametricTest = 0;
end

subplot(231); % Should match with Fig 2C
displayViolinPlot(logPowerMeditationGammaRaw,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);

subplot(232);
displayViolinPlot(logPowerMeditationGamma,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);

subplot(233); % Should match with Fig 2C
displayViolinPlot(logPowerStimulusGamma,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(234);
for pos=1:2 % 1 for meditators, 2 for controls
    data1 = logPowerMeditationGamma{pos};
    subjects1 = goodSubjectNamesMeditationGamma{pos};
    data2 = logPowerStimulusGamma{pos};
    subjects2 = goodSubjectNamesStimulusGamma{pos};
    getCorrelation(data1,subjects1,data2,subjects2);
    hold on;
end

subplot(235);
for pos=1:2 % 1 for meditators, 2 for controls
    data1 = logPowerMeditationGammaM2{pos};
    subjects1 = goodSubjectNamesMeditationGammaM2{pos};
    data2 = logPowerStimulusGammaM2{pos};
    subjects2 = goodSubjectNamesStimulusGammaM2{pos};
    getCorrelation(data1,subjects1,data2,subjects2);
    hold on;
end

% subplot(235)
% 
% practiceHoursData = load('BK1PracticeHours.mat');
% 
% data1 = cell2mat(practiceHoursData.practiceHours(2:end,2));
% subjects1 = practiceHoursData.practiceHours(2:end,1);
% data2 = logPowerMeditationGammaRaw{1};
% subjects2 = goodSubjectNamesMeditationGammaRaw{1};
% getCorrelation(data1,subjects1,data2,subjects2);
% 
% subplot(236)
% [subjects1,~,~,data1] = getDemographicDetails('BK1');
% 
% data2 = logPowerStimulusGamma{1};
% subjects2 = goodSubjectNamesStimulusGamma{1};
% getCorrelation(data1,subjects1,data2,subjects2);