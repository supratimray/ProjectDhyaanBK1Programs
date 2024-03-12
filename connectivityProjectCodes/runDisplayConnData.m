% plotFigures
comparisonStr = 'paired';
protocolName = 'G2';
analysisChoice = 'st';

badEyeCondition = 'ep';
badTrialVersion = 'v8';

freqRangeList{1} = [8 13];
freqRangeList{2} = [24 34];
freqRangeList{3} = [30 80];

axisRangeList{1} = [0 1]; axisRangeName{1} = 'YLims';
axisRangeList{2} = [0 1]; axisRangeName{2} = 'cLims (topo)';

cutoffList = [2 30];
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

electrodeList = [19 49 52]; % [14 44 47], [19 49 52], [16 17 18 48]
connMethod = 'ppc';

[connDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] = displayConnDataAllSubjects(subjectNameLists,protocolName,analysisChoice,electrodeList,connMethod,badEyeCondition,badTrialVersion,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,1);
