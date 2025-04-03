% displayPowerDataAllSubjectsMedSegmented
%
% This function analyzes and visualizes EEG power spectra across different
% meditation segments for groups of subjects. It compares between
% meditators and controls, with flexible options for reference conditions,
% frequency bands of interest, and electrode grouping.
%
% Usage:
%   [psdData, powerData, goodSubjectLists, topoplotData, freqVals] = ...
%   displayPowerDataAllSubjectsMedSegmented(hAllPlots, subjectNameLists, segmentList, ...)
%
% Parameters:
%   hAllPlots           - Plot handles structure or [] to create new plots
%   subjectNameLists    - Cell array of subject name lists {meditators, controls}
%   segmentList         - Cell a    rray of segment names to analyze
%   protocolStr         - Protocol string ('M1' or 'M2'), default 'G1'
%   analysisChoice      - Analysis type: 'st' (stimulus), 'bl' (baseline), 'combined'
%   refChoice           - Reference condition: 'none' or protocol name (e.g., 'G1')
%   badEyeCondition     - Bad eye condition type (default: 'ep')
%   badTrialVersion     - Bad trial version (default: 'v8')
%   badElectrodeRejectionFlag - How to handle bad electrodes:
%                              1: Reject bad electrodes of current protocol
%                              2: Reject common bad electrodes across all protocols
%                              3: Reject bad electrodes of G1
%   stRange             - Stimulus time range [start end] in seconds
%   freqRangeList       - Cell array of frequency ranges to analyze
%   useMedianFlag       - Use median (1) or mean (0) for averaging
%   axisRangeList       - Cell array of axis limits for plots
%   cutoffList          - Thresholds: [min electrodes, min trials, max duration]
%   pairedDataFlag      - Whether data is paired (1) or unpaired (0)
%   elecChoice          - Electrode grouping ('Occipital' or 'Frontal')
%
% Returns:
%   psdDataToReturn       - PSD data for all conditions
%   powerDataToReturn     - Power data across frequency bands
%   goodSubjectNameListsToReturn - Lists of subjects passing quality criteria
%   topoplotDataToReturn  - Topographic data for visualization
%   freqVals              - Frequency values (Hz)
%
% Related functions:
%   runDisplayMeditationSegments - GUI for meditation segment analysis
%
% Notes:
%   - This function is adapted from displayPowerDataAllSubjects
%   - The function handles paired and unpaired data differently to ensure
%     appropriate statistical comparisons

function [psdDataToReturn,powerDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] =  displayPowerDataAllSubjectsMedSegmented(hAllPlots, subjectNameLists, segmentList, protocolStr, ...
    analysisChoice, refChoice, badEyeCondition, badTrialVersion, badElectrodeRejectionFlag, ...
    stRange,freqRangeList, useMedianFlag, axisRangeList, cutoffList, pairedDataFlag, elecChoice)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Variable Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section checks and initializes variables if they don't exist or are empty
if ~exist('protocolStr','var') || isempty(protocolStr);                 protocolStr='G1';               end
if ~exist('analysisChoice','var') || isempty(analysisChoice);           analysisChoice='st';            end
if ~exist('refChoice','var') || isempty(refChoice);                     refChoice='none';               end
if ~exist('badEyeCondition','var') || isempty(badEyeCondition);         badEyeCondition='ep';           end
if ~exist('badTrialVersion','var') || isempty(badTrialVersion);         badTrialVersion='v8';           end
if ~exist('badElectrodeRejectionFlag','var');                           badElectrodeRejectionFlag=1;    end
if ~exist('stRange','var') || isempty(stRange);                         stRange = [0.25 1.25];          end
if ~exist('elecChoice','var') || isempty(elecChoice);                   elecChoice = 'Occipital';       end
if ~exist('useMedianFlag','var');                                       useMedianFlag = 0;              end
if ~exist('hAllPlots','var');                                           hAllPlots = [];                 end
if ~exist('pairedDataFlag','var');                                      pairedDataFlag = 0;             end
if ~exist('displayDataFlag','var');                                     displayDataFlag = 1;            end

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

% Get the cutoff thresholds from the cutoffList parameter
cutoffNumElectrodes = cutoffList(1);  % Minimum number of good electrodes required per subject
cutoffNumTrials     = cutoffList(2);  % Minimum number of trials required per condition
timeCutoff          = cutoffList(3);  % Maximum protocol duration for M1/M2 (in minutes) for inclusion

% Determine which electrode group to analyze based on elecChoice
if strcmp(elecChoice,'Occipital')
    % Use occipital electrode group (index 1)
    elecGroup = 1;
else
    % Use frontal electrode group (index 2) as default
    elecGroup = 2;
end

numFreqRanges = length(freqRangeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge       = 10;
displaySettings.tickLengthMedium    = [0.01 0];

% default colors
% Define color gradients for visualization

% Meditators: Red gradient from light to dark
% displaySettingsGroupMed.colorNames = [
%     0.8 0.0 0.0;  % Light red
%     0.5 0.0 0.0;  % Medium red
%     0.2 0.0 0.0   % Dark red
%     ];
% 
% % Controls: Blue gradient from light to dark
% displaySettingsGroupCon.colorNames = [
%     0.0 0.0 0.8;  % Light blue
%     0.0 0.0 0.5;  % Medium blue
%     0.0 0.0 0.2   % Dark blue
%     ];

% other color options (pink and cyan shades)
displaySettingsGroupMed.colorNames = [
    0.8026  0.3880  0.6052
    0.8732  0.1307  0.4732
    0.6039  0.1000  0.3216
];

displaySettingsGroupCon.colorNames = [
    0.0000  0.5650  0.8307
    0.0301  0.3610  0.8000
    0.0800  0.1510  0.7800
];

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
numSegments = length(segmentList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if displayDataFlag
    if isempty(hAllPlots)
        hPSD   = getPlotHandles(1,numSegments+2,[0.05 0.55 0.6 0.3],0.02,0.02,1);
        hPower = getPlotHandles(numFreqRanges, numSegments+2, [0.05 0.05 0.6 0.45], 0.02, 0.02, 0);
        hTopo  = getPlotHandles(numFreqRanges, 2*numSegments, [0.675 0.05 0.3 0.45], 0.02, 0.02, 1);
    else
        hPSD   = hAllPlots.hPSD;
        hPower = hAllPlots.hPower;
        hTopo  = hAllPlots.hTopo;
    end
    try
        x = load([capType '.mat']);
        montageChanlocs = x.chanlocs;
    catch ME
        error('Error loading channel locations: %s', ME.message);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

protocolNameList    = [{'M1a'} {'M1b'} {'M1c'} {'M2a'} {'M2b'} {'M2c'}]; % also called the "Main" protocol list
refProtocolNameList = [{'EO1'} {'EC1'} {'G1'}  {'M1'}  {'G2'}  {'EO2'} {'EC2'} {'M2'}]; % Reference protocol list containing the original 8 protocols

% Initialize arrays to collect PSD data for all segments
allSegmentsPSD = cell(2, numSegments); % {meditators/controls, segments}
goodSubjectNameListsTMP = cell(numSegments,2);

if contains(protocolStr,'M1')
    baseProtocolIndex=0;
    protocolToCheckTimeIndex = 1; % setting protocol index for using the timeCuttoff
else
    baseProtocolIndex=3;
    protocolToCheckTimeIndex = 2; % for M2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% get good subject's segmentWise %%%%%%%%%%%%%%%%

% Pre-calculate reference protocol position and flags outside the loop
isRefFromMainProtocol = 0;
refSubsegmentIndex = 0;
compareRefFlag = 0;
if ~strcmp(refChoice,'none')
    % First check if reference is in 'reference' ProtocolNameList
    protocolPosRef = find(strcmp(refProtocolNameList,refChoice));
    if isempty(protocolPosRef)
        % If not found, check in 'Main' protocolNameList (segments of M1/M2)
        protocolPosRef = find(strcmp(protocolNameList,refChoice));
        % Determine which subsegment (a/b/c) is being used as reference
        refSubsegmentIndex = mod(protocolPosRef-1,3) + 1;
        % update the flag
        isRefFromMainProtocol = 1;
        compareRefFlag = 1;
    end
else
    protocolPosRef = [];
end

% Loop through segments and get good subjects
for s=1:numSegments
    protocolName = protocolNameList{baseProtocolIndex+s};
    protocolPos = find(strcmp(protocolNameList,protocolName));
    goodSubjectNameListsTMP(s,:) = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName,baseFolderName,timeCutoff,protocolToCheckTimeIndex,isRefFromMainProtocol);
end

%%%%%%%%%%%%%%%%%%% Get common good subjects across segments %%%%%%%%%%%%%%

% Get common subjects across all segments
meditatorList = goodSubjectNameListsTMP{1,1}; % Start with first segment
controlList = goodSubjectNameListsTMP{1,2};
for s = 2:numSegments
    % Keep only subjects that appear in all segments
    meditatorList = intersect(meditatorList, goodSubjectNameListsTMP{s,1});
    controlList = intersect(controlList, goodSubjectNameListsTMP{s,2});
end
disp(['Common subjects: ' num2str(length(meditatorList)) ' meditators, ' num2str(length(controlList)) ' controls']);
commonGoodSubjectLists = {meditatorList, controlList};

%%%%%%%%%%%%%%%%%%%% Get and plot data for the segments  %%%%%%%%%%%%%%%%%%

lineNoiseRange = [48 52];
allSegmentsPower = {cell(length(segmentList), numFreqRanges), cell(length(segmentList), numFreqRanges)}; % {meditators, controls}
for s=1:numSegments
    protocolName = protocolNameList{baseProtocolIndex+s};
    protocolPos  = find(strcmp(protocolNameList,protocolName));

    displaySettings.colorNames(1,:)=displaySettingsGroupMed.colorNames(s,:);
    displaySettings.colorNames(2,:)=displaySettingsGroupCon.colorNames(s,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(commonGoodSubjectLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,saveFolderName,baseFolderName,isRefFromMainProtocol);

    %%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%
    badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
    freqPosList = cell(1,numFreqRanges); % Preallocate cell array
    for i = 1:numFreqRanges
        freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2))),badFreqPos);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numElectrodes = size(powerData{1},1);
    comparisonData = zeros(numFreqRanges,2,numElectrodes);
    topoplotDataToReturn = cell(2,numFreqRanges);

    for i=1:2 % Loop through groups (1=Meditators, 2=Controls)
        for j=1:numFreqRanges
            % Calculate power for frequency band
            freqBandPower = sum(powerData{i}(:,freqPosList{j},:),2);

            % Calculate log power or relative change
            if isempty(protocolPosRef)
                powerValues = log10(freqBandPower);
            else
                refPower = sum(powerDataRef{i}(:,freqPosList{j},:),2);
                powerValues = 10*(log10(freqBandPower) - log10(refPower));
            end

            % Get average/median power across trials
            if useMedianFlag
                data = squeeze(median(powerValues,3,'omitnan'));
            else
                data = squeeze(mean(powerValues,3,'omitnan'));
            end

            % Store results
            comparisonData(j,i,:) = data;
            topoplotDataToReturn{i,j} = data;

            % Create topoplot if display is enabled
            if displayDataFlag
                if i == 1 % Meditators
                    plotCol = s;
                else % Controls
                    plotCol = length(segmentList) + s;
                end
                set(gcf, 'CurrentAxes', hTopo(j,plotCol));
                topoplot(data, montageChanlocs, 'electrodes', 'on', ...
                    'maplimits', cLimsTopo, 'plotrad', 0.6, 'headrad', 0.6);

                % Add title for first row only
                if j==1
                    groupLabel = {'M','C'};
                    title([protocolName '-' groupLabel{i}]);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%
    psdDataToReturn = cell(1,numGroups);
    powerDataToReturn = cell(numGroups,numFreqRanges);
    goodSubjectNameListsToReturn = cell(numGroups,2);

    %%%%%%%%%%%%% Find bad subjects based on electrode cutoff %%%%%%%%%
    for i=elecGroup % numGroups
        badSubjectPosList = cell(1,2);
        for j=1:2 % meditators and controls
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

        % Get power data
        meanPSDData     = cell(1,2);
        logPSDData      = cell(1,2);
        meanPSDDataRef  = cell(1,2);

        for j=1:2 % meditators and controls
            pData = powerData{j}(electrodeGroupList{i},:,:);
            if ~isempty(protocolPosRef)
                pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
            end
            badSubjectPos = badSubjectPosList{j};

            tmp = commonGoodSubjectLists{j};
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
                set(gcf, 'CurrentAxes', hPSD(s));
                yline(0,'--k','LineWidth',2);
            end
        end

        % Violin plots for power
        for j=1:numFreqRanges
            tmpLogPower = cell(1,2);

            % Process each group (meditators and controls)
            for k=1:2
                % Calculate band power - sum across frequencies in the band
                bandPower = squeeze(sum(meanPSDData{k}(:,freqPosList{j}),2));

                if isempty(protocolPosRef)
                    % Absolute log power when no reference protocol
                    tmpLogPower{k} = log10(bandPower);
                else
                    % Calculate relative power using reference protocol
                    refBandPower = squeeze(sum(meanPSDDataRef{k}(:,freqPosList{j}),2));
                    refBandPower(refBandPower <= 0) = eps;

                    % Calculate relative change (10× difference in log scale)
                    tmpLogPower{k} = 10 * (log10(bandPower) - log10(refBandPower));
                end

                % Store data in the collection array for this segment and frequency range
                allSegmentsPower{k}{s,j} = tmpLogPower{k};
            end

            % Store band power data for return
            powerDataToReturn{i,j} = tmpLogPower;

            % Create visualization if display is enabled
            if displayDataFlag
                % Set up and display violin plot
                displaySettings.plotAxes = hPower(j,s);
                displayViolinPlot(tmpLogPower, [{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}], ...
                    1, 1, 1, pairedDataFlag, displaySettings);

                % Add appropriate labels
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

%%%%%%%%%%%%%%%%%  plot power data for all the segments together %%%%%%%%%%
if displayDataFlag
    % Define color maps for Meditators and Controls
    colorMapMed = num2cell(displaySettingsGroupMed.colorNames, 2);
    colorMapCon = num2cell(displaySettingsGroupCon.colorNames, 2);

    % Plot violin plots for each frequency range, looping through groups
    for j=1:numFreqRanges
        for g=1:2
            if g==1
                groupData = allSegmentsPower{1}(:,j)';
                colorMap  = colorMapMed;
                axisIndex = 4;
                groupName = 'Meditators';
            else
                groupData = allSegmentsPower{2}(:,j)';
                colorMap  = colorMapCon;
                axisIndex = 5;
                groupName = 'Controls';
            end

            displaySettings.plotAxes = hPower(j,axisIndex);
            displayViolinPlot(groupData, colorMap, 1, 1, 1, 1, displaySettings);
            if j==1, title(hPower(j,axisIndex),groupName); end
        end
    end

    % Prepare arrays for PSD plots of meditators and controls
    groupData      = {allSegmentsPSD(1,:), allSegmentsPSD(2,:)};
    groupColors    = {cat(1, colorMapMed{:}), cat(1, colorMapCon{:})};
    groupTitles    = {'Meditators - All Segments','Controls - All Segments'};
    groupAxes      = [hPSD(4) hPSD(5)];

    % Loop through groups for PSD plots
    for g=1:2
        set(gcf, 'CurrentAxes', groupAxes(g));
        displaySettingsGroup       = displaySettings;
        displaySettingsGroup.colorNames = groupColors{g};

        displayAndcompareData(groupAxes(g), groupData{g}, freqVals, ...
            displaySettingsGroup, yLimsPSD, 1, useMedianFlag, ...
            ~pairedDataFlag, compareRefFlag, refSubsegmentIndex);

        title(groupTitles{g});
        xlim(freqLims); ylim(yLimsPSD);
        if ~isempty(protocolPosRef), yline(0,'--k','LineWidth',2); end
    end
end
end  % main function end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            disp(['Timing cutoff exceeded for subject: ' subjectName]);
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

%%%%%%%%%%%%%%%%%%%%%%%%% getPowerDataAllSubjects %%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%% getPowerData %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%% getBadElectrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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