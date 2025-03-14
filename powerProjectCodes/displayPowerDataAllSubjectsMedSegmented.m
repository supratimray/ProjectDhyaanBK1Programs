% analysisChoice - 'st', 'bl' or 'combined'
% refChoice - 'none' (show raw PSDs and power) or a protocolName. For
% example, if refChoice is 'G1' then we use G1 baseline as reference.

% badTrialRejectionFlag:
% 1: Reject badElectrodes of protocolName
% 2. Reject common badElectrodes of all protocols
% 3: Reject badElectrodes of G1

% Option added to return PSD, power and topoplot data. Also to simply
% return these without displaying here

% adapted and modified from 'displayPowerDataAllSubjects'

function [psdDataToReturn,powerDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] =  displayPowerDataAllSubjectsMedSegmented(hAllPlots, subjectNameLists, segmentList, ~, protocolStr, ...
    analysisChoice, refChoice, badEyeCondition, badTrialVersion, badElectrodeRejectionFlag, ...
    freqRangeList, useMedianFlag, axisRangeList, cutoffList, pairedDataFlag, elecChoice)

if ~exist('protocolStr','var');                 protocolStr='G1';               end
if ~exist('analysisChoice','var');              analysisChoice='st';            end
if ~exist('refChoice','var');                   refChoice='none';               end
if ~exist('badEyeCondition','var');             badEyeCondition='ep';           end
if ~exist('badTrialVersion','var');             badTrialVersion='v8';           end
if ~exist('badElectrodeRejectionFlag','var');   badElectrodeRejectionFlag=1;    end
if ~exist('stRange','var');                     stRange = [0.25 1.25];          end
if ~exist('elecChoice','var');                  elecChoice = 'Occipital';       end
if ~exist('useMedianFlag','var');               useMedianFlag = 0;              end
if ~exist('hAllPlots','var');                   hAllPlots = [];                 end
if ~exist('pairedDataFlag','var');              pairedDataFlag = 0;             end
if ~exist('displayDataFlag','var');             displayDataFlag = 1;            end

if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13];  % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [3 30 20];
end
cutoffNumElectrodes = cutoffList(1);
cutoffNumTrials     = cutoffList(2);
timeCuttOff         = cutoffList(3);

if strcmp(elecChoice,'Occipital')
    elecGroup           = 1;
else
    elecGroup           = 2; % frontal
end

numFreqRanges = length(freqRangeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge       = 10;
displaySettings.tickLengthMedium    = [0.01 0];

% default colors
displaySettings.colorNames(1,:)     = [0.8 0 0.8];      % Purple
displaySettings.colorNames(2,:)     = [0.25 0.41 0.88]; % Cyan
displaySettingsGroupMed.colorNames  = [0.8  0  0   ; 0.5 0  0   ;  0.2  0  0  ];     % Light to dark red shades
displaySettingsGroupCon.colorNames  = [0    0  0.8 ; 0   0  0.5 ;  0    0  0.2];

titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

freqLims = axisRangeList{1};
yLimsPSD = axisRangeList{2};
cLimsTopo = axisRangeList{3};

%%%%%%%%%%%%%%%%%%%%%%%% Get electrode groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType  = 'actiCap64_UOL';
baseFolderName = 'savedData';
saveFolderName = 'savedDataMedSegmented';

[electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
numGroups = length(electrodeGroupList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if displayDataFlag
    if isempty(hAllPlots)
        hPSD   = getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.3],0.02,0.02,1);
        hPower = getPlotHandles(numFreqRanges, length(segmentList)*2, [0.05 0.05 0.6 0.45], 0.02, 0.02, 0);
        hTopo  = getPlotHandles(numFreqRanges, length(segmentList)*2, [0.675 0.05 0.3 0.45], 0.02, 0.02, 1);
    else
        hPSD   = hAllPlots.hPSD;
        hPower = hAllPlots.hPower;
        hTopo  = hAllPlots.hTopo;
    end
    x = load([capType '.mat']); montageChanlocs = x.chanlocs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList    = [{'M1a'} {'M1b'} {'M1c'} {'M2a'} {'M2b'} {'M2c'}];
refProtocolNameList = [{'EO1'} {'EC1'} {'G1'}  {'M1'}  {'G2'}  {'EO2'} {'EC2'} {'M2'}];

% Initialize arrays to collect PSD data for all segments
allSegmentsPSD = cell(2, length(segmentList)); % {meditators/controls, segments}
goodSubjectNameListsTMP = cell(length(segmentList),2);

if contains(protocolStr,'M1')
    protocolBaseIndex=0;
    protocolToCheckIndex = 1; % check timeCuttoff
else
    protocolBaseIndex=3;
    protocolToCheckIndex = 2; % for M2
end

% get good subject's segmentWise
for s=1:length(segmentList)
    protocolName = protocolNameList{protocolBaseIndex+s};
    protocolPos = find(strcmp(protocolNameList,protocolName));
    if ~strcmp(refChoice,'none')
        % First check if reference is in refProtocolNameList
        refProtocolPosRefIndex = find(strcmp(refProtocolNameList,refChoice));
        isRefFromMainProtocol = 0;
        refSubsegmentIndex = 0;
        if isempty(refProtocolPosRefIndex)
            % If not found, check in protocolNameList
            refProtocolPosRefIndex = find(strcmp(protocolNameList,refChoice));
            isRefFromMainProtocol = 1;
            compareRefFlag = 1;
            % Determine which subsegment (a/b/c) is being used as reference
            refSubsegmentIndex = mod(refProtocolPosRefIndex-1,3) + 1;
        end
        protocolPosRef = refProtocolPosRefIndex;
    else
        protocolPosRef = [];
        isRefFromMainProtocol = 0;
        refSubsegmentIndex = 0;
        compareRefFlag = 0;
    end
    goodSubjectNameListsTMP(s,:) = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName,baseFolderName,timeCuttOff,protocolToCheckIndex,isRefFromMainProtocol);
end

%%%%%%%%%%%%%%%%%%% Get common subjects across segments %%%%%%%%%%%%%%%%%%%%%
% start with first segment's subjects
meditatorList = goodSubjectNameListsTMP{1,1};
controlList   = goodSubjectNameListsTMP{1,2};
for s=2:size(goodSubjectNameListsTMP,1)
    % find common subject's across all segments
    meditatorList = intersect(meditatorList, goodSubjectNameListsTMP{s,1});
    controlList   = intersect(controlList, goodSubjectNameListsTMP{s,2});
end
goodSubjectNameListsGrouped = {meditatorList, controlList};

% Initialize data collection arrays before segment loop
allSegmentsPower = {cell(length(segmentList), numFreqRanges), cell(length(segmentList), numFreqRanges)}; % {meditators, controls}

%%%%%%%%%%%%%%%%%%%% Get and plot data for the segments  %%%%%%%%%%%%%%%%%%%%%
for s=1:length(segmentList)
    protocolName = protocolNameList{protocolBaseIndex+s};
    protocolPos  = find(strcmp(protocolNameList,protocolName));

    displaySettings.colorNames(1,:)=displaySettingsGroupMed.colorNames(s,:);
    displaySettings.colorNames(2,:)=displaySettingsGroupCon.colorNames(s,:);

    if ~strcmp(refChoice,'none')
        % First check if reference is in refProtocolNameList
        refProtocolPosRefIndex = find(strcmp(refProtocolNameList,refChoice));
        isRefFromMainProtocol = 0;
        if isempty(refProtocolPosRefIndex)
            % If not found, check in protocolNameList
            isRefFromMainProtocol = 1;
            compareRefFlag = 1;
            refProtocolPosRefIndex = find(strcmp(protocolNameList,refChoice));
        end
        protocolPosRef = refProtocolPosRefIndex;
    else
        protocolPosRef = [];
        isRefFromMainProtocol = 0;
        compareRefFlag = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    goodSubjectNameLists = goodSubjectNameListsGrouped;
    [powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(goodSubjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,saveFolderName,baseFolderName,isRefFromMainProtocol);

    %%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%
    freqPosList = cell(1,numFreqRanges);
    lineNoiseRange = [48 52];
    badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
    for i = 1:numFreqRanges
        freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2))),badFreqPos);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Show Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numElectrodes = size(powerData{1},1);
    comparisonData = zeros(numFreqRanges,2,numElectrodes);
    topoplotDataToReturn = cell(2,numFreqRanges);

    for i=1:2
        for j=1:numFreqRanges
            if isempty(protocolPosRef)
                x = log10(sum(powerData{i}(:,freqPosList{j},:),2));
            else
                x = 10*(log10(sum(powerData{i}(:,freqPosList{j},:),2)) - log10(sum(powerDataRef{i}(:,freqPosList{j},:),2)));
            end
            if useMedianFlag
                data = squeeze(median(x,3,'omitnan'));
            else
                data = squeeze(mean(x,3,'omitnan'));
            end
            comparisonData(j,i,:) = data;
            topoplotDataToReturn{i,j} = data;
            if displayDataFlag
                if i == 1 % Meditators
                    plotCol = s;
                else % Controls
                    plotCol = length(segmentList) + s;
                end
                axes(hTopo(j,plotCol));
                topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo,'plotrad',0.6,'headrad',0.6);
                if j==1        % Add titles only for first row
                    if i == 1
                        title([protocolName '-M']);
                    else
                        title([protocolName '-C']);
                    end
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    psdDataToReturn = cell(1,numGroups);
    powerDataToReturn = cell(numGroups,numFreqRanges);
    goodSubjectNameListsToReturn = cell(numGroups,2);

    for i=elecGroup % numGroups
        %%%%%%%%%%%%% Find bad subjects based on electrode cutoff %%%%%%%%%%%%%
        badSubjectPosList = cell(1,2);
        for j=1:2
            pData = powerData{j}(electrodeGroupList{i},:,:);
            numGoodElecs = length(electrodeGroupList{i}) - sum(isnan(squeeze(pData(:,1,:))),1);
            badSubjectPosList{j} = find(numGoodElecs<cutoffNumElectrodes);

            if ~isempty(protocolPosRef)
                pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
                numGoodElecsRef = length(electrodeGroupList{i}) - sum(isnan(squeeze(pDataRef(:,1,:))),1);
                badSubjectPosRef = find(numGoodElecsRef<=cutoffNumElectrodes);
                badSubjectPosList{j} = unique(cat(2,badSubjectPosList{j},badSubjectPosRef));
            end
        end

        % if paired dataset, take common bad subjects
        if pairedDataFlag
            badSubjectPosList{1} = union(badSubjectPosList{1},badSubjectPosList{2});
            badSubjectPosList{2} = badSubjectPosList{1};
        end

        % Get data
        meanPSDData = cell(1,2);
        meanPSDDataRef = cell(1,2);
        logPSDData = cell(1,2);

        for j=1:2
            pData = powerData{j}(electrodeGroupList{i},:,:);
            if ~isempty(protocolPosRef)
                pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
            end
            badSubjectPos = badSubjectPosList{j};

            tmp = goodSubjectNameLists{j};
            tmp(badSubjectPos) = [];
            goodSubjectNameListsToReturn{i,j} = tmp;

            if ~isempty(badSubjectPos)
                disp([groupNameList{i} ', ' titleStr{j} ', '  'Not enough good electrodes for ' num2str(length(badSubjectPos)) ' subjects.']);
                pData(:,:,badSubjectPos)=[];
                if ~isempty(protocolPosRef)
                    pDataRef(:,:,badSubjectPos)=[];
                end
            end
            meanPSDData{j} = squeeze(mean(pData,1,'omitnan'))';

            if isempty(protocolPosRef)
                logPSDData{j} = log10(meanPSDData{j});
            else
                meanPSDDataRef{j} = squeeze(mean(pDataRef,1,'omitnan'))';
                logPSDData{j} = 10*(log10(meanPSDData{j}) - log10(meanPSDDataRef{j}));
            end

            if displayDataFlag
                text(30,yLimsPSD(2)-0.5*j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(s));
            end
        end

        psdDataToReturn{i} = logPSDData;

        % Store PSD data for each segment
        allSegmentsPSD{1,s} = logPSDData{1}; % Meditators
        allSegmentsPSD{2,s} = logPSDData{2}; % Controls
        if displayDataFlag
            displayAndcompareData(hPSD(s),logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,~pairedDataFlag);
            ylim(hPSD(s),yLimsPSD);
            % axis tight
            title(hPSD(s),[groupNameList{i} ' - Segment ' segmentList{s}]);  % Added segment name
            xlim(hPSD(s),freqLims);

            if ~isempty(protocolPosRef)
                axes(hPSD(s));
                yline(0,'--k','LineWidth',2);
            end
        end

        % Violin plots for power
        for j=1:numFreqRanges
            tmpLogPower = cell(1,2);
            for k=1:2
                if isempty(protocolPosRef)
                    tmpLogPower{k} = log10(squeeze(sum(meanPSDData{k}(:,freqPosList{j}),2)));
                else
                    tmpLogPower{k} = 10*(log10(squeeze(sum(meanPSDData{k}(:,freqPosList{j}),2))) - log10(squeeze(sum(meanPSDDataRef{k}(:,freqPosList{j}),2))));
                end

                % Store data in the collection array
                allSegmentsPower{k}{s,j} = tmpLogPower{k};
            end
            powerDataToReturn{i,j} = tmpLogPower;
            if displayDataFlag
                % display violin plots for power
                displaySettings.plotAxes = hPower(j,s);
                displayViolinPlot(tmpLogPower,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);
                if i==1
                    % Add frequency info to y-label
                    ylabel(hPower(j,i),['\Delta Power (' num2str(freqRangeList{j}(1)) '-' num2str(freqRangeList{j}(2)) ' Hz)']);
                elseif j==1
                    title(hPower(j,s),[groupNameList{i} ' - Segment ' segmentList{s}]);  % Added segment name
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%  plot power data for all the segments together %%%%%%%%%%%%%%%%%
if displayDataFlag
    for j=1:numFreqRanges
        combinedDataMed = allSegmentsPower{1}(:,j)';
        combinedDataCon = allSegmentsPower{2}(:,j)';
        % Plot meditators
        displaySettings.plotAxes = hPower(j,4);
        displaySettings.xTickLabels = segmentList;
        displaySettings.parametricTest = 1;
        displayViolinPlot(combinedDataMed, [{[0.8 0 0]} {[0.5 0 0]} {[0.2 0 0]}], ...
            1, 1, 1, 1, displaySettings);
        % Plot controls
        displaySettings.plotAxes = hPower(j,5);
        displayViolinPlot(combinedDataCon, {[0 0 0.8], [0 0 0.5], [0 0 0.2]}, ...
            1, 1, 1, 1, displaySettings);
        if j==1
            title(hPower(j,4),'Meditators');
            title(hPower(j,5),'Controls');
        end
    end

    % Plot meditator's PSDs in column 4
    axes(hPSD(4));
    displaySettingsGroupMed = displaySettings;
    displaySettingsGroupMed.colorNames = [0.8  0   0;   0.5 0  0;    0.2  0    0];     % Light to dark red shades
    displayAndcompareData(hPSD(4),allSegmentsPSD(1,:),freqVals,displaySettingsGroupMed,yLimsPSD,1,useMedianFlag,~pairedDataFlag,compareRefFlag,refSubsegmentIndex);
    title('Meditators - All Segments');
    xlim(freqLims);
    ylim(yLimsPSD);
    if ~isempty(protocolPosRef)
        yline(0,'--k','LineWidth',2);
    end

    % Plot control's PSDs in column 5
    axes(hPSD(5));
    cla;
    displaySettingsGroupCon = displaySettings;
   displaySettingsGroupCon.colorNames = [0 0 0.8 ; 0 0 0.5 ; 0 0 0.2];
    displayAndcompareData(hPSD(5),allSegmentsPSD(2,:),freqVals,displaySettingsGroupCon,yLimsPSD,1,useMedianFlag,~pairedDataFlag,compareRefFlag,refSubsegmentIndex);
    title('Controls - All Segments');
    xlim(freqLims);
    ylim(yLimsPSD);
    if ~isempty(protocolPosRef)
        yline(0,'--k','LineWidth',2);
    end
end
end  % main function end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName,baseFolderName,timeCuttOff,protcolToCheckIndex,isRefFromMainProtocol)

% For unpaired case, subjects can be rejected if either data in analysis or
% ref period is bad. For paired case, a pair is rejected even if one of the two
% subjects in the pair is bad. Based on the condition, we get a new good
% subject list. For both the cases, a subject is rejected if the time exceeds the Cut-off timing

badSubjectIndex = cell(1,2);
badSubjectIndexRef = cell(1,2);

for i=1:2
    numSubjects = length(subjectNameLists{i});
    badSubjectIndexTMP = zeros(1,numSubjects);
    badSubjectIndexRefTMP = zeros(1,numSubjects);

    for j=1:numSubjects
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_medseg.mat']));

        if isempty(getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials))
            disp(['Not enough trials for subject: ' subjectName]);
            badSubjectIndexTMP(j)=1;
        end

        % use the timeCuffOff to find the badSubjects
        protocolLengthTimes = tmpData.protcolLengthTime;
        protocolLengthTime = protocolLengthTimes(protcolToCheckIndex);

        if protocolLengthTime>timeCuttOff
            badSubjectIndexTMP(j)=1;
        end

        if ~isempty(protocolPosRef)
            if isRefFromMainProtocol
                % If reference is from protocolNameList, check from same file
                if isempty(getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials))
                    disp(['Not enough trials in ref period for subject: ' subjectName]);
                    badSubjectIndexRefTMP(j)=1;
                end
            else
                % If reference is from refProtocolNameList, load from baseFolderName
                tmpData2 = load(fullfile(baseFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));
                if isempty(getPowerData(tmpData2,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials))
                    disp(['Not enough trials in ref period for subject: ' subjectName]);
                    badSubjectIndexRefTMP(j)=1;
                end
            end
        end
    end
    badSubjectIndex{i} = badSubjectIndexTMP;
    badSubjectIndexRef{i} = badSubjectIndexRefTMP;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Now find good subjects %%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectNameLists = cell(1,2);

if ~pairedDataFlag
    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        if isempty(protocolPosRef)
            badPos = find(badSubjectIndex{i});
        else
            badPos = union(find(badSubjectIndex{i}),find(badSubjectIndexRef{i}));
        end
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
else

    if isempty(protocolPosRef)
        badPos = find(sum(cell2mat(badSubjectIndex')));

    else
        badPos = union(find(sum(cell2mat(badSubjectIndex'))),find(sum(cell2mat(badSubjectIndexRef'))));
    end

    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
end
end
function [powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,saveFolderName,baseFolderName,isRefFromMainProtocol)

powerData = cell(1,2);
powerDataRef = cell(1,2);

for i=1:2
    powerDataTMP=[];
    powerDataRefTMP=[];

    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_medseg.mat'] ));
        freqVals = tmpData.freqVals;

        tmpPower = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials);
        powerDataTMP = cat(3,powerDataTMP,tmpPower);

        if ~isempty(protocolPosRef)
            if isRefFromMainProtocol
                % If reference is from protocolNameList, use the same segmented data file
                tmpPowerRef = getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
            else
                % If reference is from refProtocolNameList, load from baseFolderName
                tmpData2 = load(fullfile(baseFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));
                tmpPowerRef = getPowerData(tmpData2,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
            end
            powerDataRefTMP = cat(3,powerDataRefTMP,tmpPowerRef);
        end
    end
    powerData{i} = powerDataTMP;
    powerDataRef{i} = powerDataRefTMP;
end
end
function tmpPower = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials)

numTrials = tmpData.numTrials(protocolPos);
badElectrodes = getBadElectrodes(tmpData.badElectrodes,badElectrodeRejectionFlag,protocolPos);

if numTrials < cutoffNumTrials
    tmpPower = [];
else
    if strcmpi(analysisChoice,'st')
        tmpPower = tmpData.psdValsST{protocolPos};
    elseif strcmpi(analysisChoice,'bl')
        tmpPower = tmpData.psdValsBL{protocolPos};
    else
        tmpPower = (tmpData.psdValsST{protocolPos}+tmpData.psdValsBL{protocolPos})/2; % average
    end
    tmpPower(badElectrodes,:) = NaN;
end
end
function badElectrodes = getBadElectrodes(badElectrodeList,badElectrodeRejectionFlag,protocolPos)

if badElectrodeRejectionFlag==1 % Bad electrodes for the protocol
    badElectrodes = badElectrodeList{protocolPos};
elseif badElectrodeRejectionFlag==2 % common bad electrodes for all protocols
    badElectrodes=[];
    for i=1:length(badElectrodeList)
        badElectrodes=cat(1,badElectrodes,badElectrodeList{i});
    end
    badElectrodes = unique(badElectrodes);
elseif badElectrodeRejectionFlag==3 % Bad electrodes of G1
    badElectrodes = badElectrodeList{3};
end
end