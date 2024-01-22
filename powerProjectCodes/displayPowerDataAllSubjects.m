% analysisChoice - 'st', 'bl' or 'combined'

% refChoice - 'none' (show raw PSDs and power) or a protocolName. For
% example, if refChoice is 'G1' then we use G1 baseline as reference.

% badTrialRejectionFlag:
% 1: Reject badElectrodes of protocolName
% 2. Reject common badElectrodes of all protocols
% 3: Reject badElectrodes of G1

% Option added to return PSD, power and topoplot data. Also to simply
% return these without displaying here.

function [psdDataToReturn,powerDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] = displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,displayDataFlag)

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('refChoice','var');             refChoice='none';             end

if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=1;  end

if ~exist('stRange','var');               stRange = [0.25 1.25];        end

if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [3 30];
end
cutoffNumElectrodes = cutoffList(1);
cutoffNumTrials = cutoffList(2);

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [1 0 0];
displaySettings.colorNames(2,:) = [0 1 0];
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

freqLims = axisRangeList{1};
yLimsPSD = axisRangeList{2};
cLimsTopo = axisRangeList{3};

%%%%%%%%%%%%%%%%%%%%%%%% Get electrode groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';
saveFolderName = 'savedData';

[electrodeGroupList,groupNameList] = getElectrodeGroups(gridType,capType);
numGroups = length(electrodeGroupList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    if isempty(hAllPlots)
        hPSD  = getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.3],0.02,0.02,1);
        hPower = getPlotHandles(numFreqRanges,numGroups,[0.05 0.05 0.6 0.45],0.02,0.02,0);
        hTopo0 = getPlotHandles(1,2,[0.675 0.7 0.3 0.15],0.02,0.02,1);
        hTopo1 = getPlotHandles(1,3,[0.675 0.55 0.3 0.13],0.02,0.02,1);
        hTopo2 = getPlotHandles(numFreqRanges,3,[0.675 0.05 0.3 0.45],0.02,0.02,1);
    else
        hPSD = hAllPlots.hPSD;
        hPower = hAllPlots.hPower;
        hTopo0 = hAllPlots.hTopo0;
        hTopo1 = hAllPlots.hTopo1;
        hTopo2 = hAllPlots.hTopo2;
    end
    montageChanlocs = showElectrodeGroups(hTopo0(1,:),capType,electrodeGroupList,groupNameList);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
protocolPos = find(strcmp(protocolNameList,protocolName));

if ~strcmp(refChoice,'none')
    protocolPosRef = find(strcmp(protocolNameList,refChoice));
else
    protocolPosRef = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName);
[powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(goodSubjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,saveFolderName);

%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
lineNoiseRange = [48 52];
badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
for i = 1:numFreqRanges
    freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2))),badFreqPos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElectrodes = size(powerData{1},1);
percentData = zeros(2,numElectrodes);
comparisonData = zeros(numFreqRanges,2,numElectrodes);

topoplotDataToReturn = cell(2,numFreqRanges);
for i=1:2
    if isempty(protocolPosRef)
        x=powerData{i};
    else
        x=powerData{i} ./ powerDataRef{i};
    end
    numSubjects = size(x,3);

    %%%%%%%%%%%% Show percent of bad subjects per electrode %%%%%%%%%%%%%%%
    numBadSubjects = zeros(1,numElectrodes);
    for j=1:numElectrodes
        numBadSubjects(j) = sum(isnan(squeeze(x(j,1,:))));
    end

    if displayDataFlag
        axes(hTopo1(i)); %#ok<*LAXES>

        % Modification in the topoplot code which allows us to not interpolate across electrodes.
        %topoplot_murty(numBadSubjects/numSubjects,montageChanlocs,'electrodes','off','style','blank','drawaxis','off','emarkercolors',numBadSubjects/numSubjects); colorbar;
        percentData(i,:) = 100*(numBadSubjects/numSubjects);
        topoplot(percentData(i,:),montageChanlocs,'maplimits',[0 100],'electrodes','on','plotrad',0.6,'headrad',0.6); colorbar;
        title(titleStr{i},'color',displaySettings.colorNames(i,:));
        if i==1
            ylabel('Bad subjects (%)');
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%% Show topoplots of power %%%%%%%%%%%%%%%%%%%%%%%
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
            axes(hTopo2(j,i));
            topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo,'plotrad',0.6,'headrad',0.6); colorbar;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plot the difference of topoplots %%%%%%%%%%%%%%%%%%%
if displayDataFlag
    axes(hTopo1(3));
    topoplot(-diff(percentData),montageChanlocs,'maplimits',[-25 25],'electrodes','on','plotrad',0.6,'headrad',0.6); colorbar;

    for i=1:numFreqRanges
        axes(hTopo2(i,3));
        data = -diff(squeeze(comparisonData(i,:,:)));
        if isempty(protocolPosRef)
            data = 10*data;
        end
        topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo,'plotrad',0.6,'headrad',0.6); colorbar;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psdDataToReturn = cell(1,numGroups);
powerDataToReturn = cell(numGroups,numFreqRanges);
goodSubjectNameListsToReturn = cell(numGroups,2);

for i=1:numGroups

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
            text(30,yLimsPSD(2)-0.5*j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(i));
        end
    end

    psdDataToReturn{i} = logPSDData;
    if displayDataFlag
        displayAndcompareData(hPSD(i),logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,~pairedDataFlag);
        title(groupNameList{i});
        xlim(hPSD(i),freqLims);
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
        end

        powerDataToReturn{i,j} = tmpLogPower;

        if displayDataFlag
            % display violin plots for power
            displaySettings.plotAxes = hPower(j,i);
            if ~useMedianFlag
                displaySettings.parametricTest = 1;
            else
                displaySettings.parametricTest = 0;
            end

            if i==numGroups && j==1
                displaySettings.showYTicks=1;
                displaySettings.showXTicks=1;
            else
                displaySettings.showYTicks=0;
                displaySettings.showXTicks=0;
            end
            displayViolinPlot(tmpLogPower,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);
            if i==1
                ylabel(hPower(j,i),[num2str(freqRangeList{j}(1)) '-' num2str(freqRangeList{j}(2)) ' Hz'],'color',freqRangeColors(j,:));
            end

            % Add lines in PSD plots
            for k=1:2
                line([freqRangeList{j}(k) freqRangeList{j}(k)],yLimsPSD,'color',freqRangeColors(j,:),'parent',hPSD(i));
            end

            if ~isempty(protocolPosRef)
                line([0 freqVals(end)],[0 0],'color','k','parent',hPSD(i));
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,saveFolderName)

% For unpaired case, subjects can be rejected if either data in analysis or
% ref period is bad. For paired case, a pair is rejected even if one of the two
% subjects in the pair is bad. Based on the condition, we get a new good
% subject list.

badSubjectIndex = cell(1,2);
badSubjectIndexRef = cell(1,2);

for i=1:2
    numSubjects = length(subjectNameLists{i});
    badSubjectIndexTMP = zeros(1,numSubjects);
    badSubjectIndexRefTMP = zeros(1,numSubjects);

    for j=1:numSubjects
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));

        if isempty(getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials))
            disp(['Not enough trials for subject: ' subjectName]);
            badSubjectIndexTMP(j)=1;
        end

        if ~isempty(protocolPosRef)
            if isempty(getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials))
                disp(['Not enough trials in ref period for subject: ' subjectName]);
                badSubjectIndexRefTMP(j)=1;
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
function [powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,saveFolderName)

powerData = cell(1,2);
powerDataRef = cell(1,2);

for i=1:2
    powerDataTMP=[];
    powerDataRefTMP=[];

    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));
        freqVals = tmpData.freqVals;

        tmpPower = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials);
        powerDataTMP = cat(3,powerDataTMP,tmpPower);

        if ~isempty(protocolPosRef)
            tmpPowerRef = getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
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
function montageChanlocs = showElectrodeGroups(hPlots,capType,electrodeGroupList,groupNameList)

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']); montageChanlocs = x.chanlocs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hPlots(1));
electrodeSize = 5;
numElectrodes = length(montageLabels);

for i=1:numElectrodes
    montageChanlocs(i).labels = ' ';
end

numElectrodeGroups = length(electrodeGroupList);
electrodeGroupColorList = jet(numElectrodeGroups);

for i=1:numElectrodeGroups
    topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize},'plotrad',0.6,'headrad',0.6);
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank','plotrad',0.6,'headrad',0.6);

axes(hPlots(2))
set(hPlots(2),'visible','off');
for i=1:numElectrodeGroups
    text(0.05,0.9-0.15*(i-1),[groupNameList{i} '(N=' num2str(length(electrodeGroupList{i})) ')'],'color',electrodeGroupColorList(i,:),'unit','normalized');
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
