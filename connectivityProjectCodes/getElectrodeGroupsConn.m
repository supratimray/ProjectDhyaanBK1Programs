% Returns the group of electrodes used for averaging connectivity. When
% groupType is set to 'abs', a fixed set of electrode groups is returned.
% When groupType is set to 'rel', electrode groups are made relative to the
% position of the refElectrodes. A separate set of electrode groups is made
% for each reference group.

function [electrodeGroupList,groupNameList,binnedCenters,montageChanlocs] = getElectrodeGroupsConn(groupType,refElectrodes,capType)

x = load([capType '.mat']); 
montageChanlocs = x.chanlocs;

if strcmp(groupType,'abs')
    [electrodeGroupList,groupNameList] = getElectrodeGroups('EEG',capType,0,'conn');
    binnedCenters = [];
else
    [electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConnRelative(refElectrodes,montageChanlocs);
end
end

function [electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConnRelative(refElectrodes,montageChanlocs)

%%%%% Discretize connectivty into bins depending on distance from seed %%%%
binWidth = 0.25;
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binnedCenters = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(montageChanlocs);

numElectrodes = length(refElectrodes);
electrodeGroupList = cell(numElectrodes,nbins);
groupNameList = cell(1,nbins);

for e=1:numElectrodes
    dist = sqrt((angl_dist(loc.azi(refElectrodes(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(refElectrodes(e)),loc.ele,'e')).^2);
    fitx = cos((dist/180)*pi);
    binned_fitx = discretize(fitx,binEdges);
    
    for b = 1:nbins % number of bins
        electrodeGroupList{e,nbins-b+1} = find(binned_fitx == b);
        if e==1
            groupNameList{nbins-b+1} = [num2str(binEdges(b)) '<x<' num2str(binEdges(b+1))];
        end
    end
end
end

%%%%%%%% Codes written by Santosh for TLSA connectivity project %%%%%%%%%%%
function loc = getElecLocAngles(chanlocs)
azi = zeros(1,length(chanlocs)); ele = zeros(1,length(chanlocs));
for e = 1:length(chanlocs)
    azi(e) = chanlocs(e).sph_theta;
    ele(e) = chanlocs(e).sph_phi;
end
loc.azi = azi;
loc.ele = ele;
end
function out_theta = angl_dist(in_theta_ref,in_theta,val)
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    if(in_theta_ref > 90)
    in_theta(14) = 90;
    elseif(in_theta_ref < -90)
        in_theta(14) = -90;
    end
end
in_theta = abs(in_theta-in_theta_ref);
out_theta = zeros(1,length(in_theta));
for i=1:length(in_theta)
    if(in_theta(i) > 180)
        out_theta(i) = 360 - in_theta(i); % to get shortest angular distance
    else
        out_theta(i) = in_theta(i);
    end
end
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    out_theta(14) = 0;
end
end