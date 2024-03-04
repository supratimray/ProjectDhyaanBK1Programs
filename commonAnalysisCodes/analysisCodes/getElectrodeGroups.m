function [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType,useFocalOccipitalFlag)

if ~exist('useFocalOccipitalFlag','var');    useFocalOccipitalFlag=0;        end

[~,~,~,~,~,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Occipital Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 1: Choose highPriority electrodes
if useFocalOccipitalFlag
    electrodeGroupList{1} =[13  32+[12 20 16 15 17]];
else
    electrodeGroupList{1} = highPriorityElectrodeNums; % now called occipital
end
groupNameList{1} = 'Occipital';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fronto-Temporal Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 3: Choose frontal, fronto-central and fronto-temporal electrodes
subFrontal            = [1 32      32+[1 2 29 30 5 28]];
subFrontoCental       = [6 8 25 28 32+[9 24]];
subFrontoTemporal     = [4 31 9 26 32+[6 27 10 23]];
electrodeGroupList{2} = [subFrontal subFrontoCental subFrontoTemporal]; % Fronto-Central and Frontal
groupNameList{2}      = 'Fronto-Temporal';

end