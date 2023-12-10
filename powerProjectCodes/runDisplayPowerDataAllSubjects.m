% runDisplayPowerDataAllSubjects

[~, meditatorList, controlList] = getGoodSubjectsBK1;

subjectNameLists{1} = meditatorList;
subjectNameLists{2} = controlList;

protocolName = 'G1'; analysisChoice = 'st';
displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice);