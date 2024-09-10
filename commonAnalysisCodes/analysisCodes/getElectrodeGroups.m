function [electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType,useFocalOccipitalFlag,analysisType)

if ~exist('useFocalOccipitalFlag','var');   useFocalOccipitalFlag=0;    end
if ~exist('analysisType','var');            analysisType='power';       end

if strcmpi(analysisType,'power')

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

elseif strcmpi(analysisType,'conn')
    electrodeGroupList{1} = [14:16 32+ (12:15)];                           groupNameList{1} = 'Occipital_L'; % Occipital
    electrodeGroupList{2} = [18:20 32+ (17:20)];                           groupNameList{2} = 'Occipital_R'; % Occipital

    electrodeGroupList{3} = [6:8 11 12  32+[7:9 11]];                      groupNameList{3} = 'Central_L'; % Central
    electrodeGroupList{4} = [22 23 25 28 29 32+[22 24:26]];                groupNameList{4} = 'Central_R'; % Central

%    electrodeGroupList{5} = [5 9 10 32+[6 10]];                            groupNameList{5} = 'Temporal_L'; % Temporal
%    electrodeGroupList{6} = [21 26 27 32+[23 27]];                         groupNameList{6} = 'Temporal_R'; % Temporal

    electrodeGroupList{5} = [1 3 4 32+[1 2 4 5]];                          groupNameList{5} = 'Frontal_L'; % Frontal
    electrodeGroupList{6} = [30:32 32+[28:31]];                            groupNameList{6} = 'Frontal_R'; % Frontal
end

end