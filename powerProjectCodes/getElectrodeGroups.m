function [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType)

[~,~,~,electrodeGroupList0,~,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);

% Combine some groups
electrodeGroupList{1} = highPriorityElectrodeNums; % now called occipital
groupNameList{1} = 'Occipital';
electrodeGroupList{2} = [electrodeGroupList0{3} electrodeGroupList0{4}]; % Fronto-Central and Frontal
groupNameList{2} = 'Frontal-Central';

end