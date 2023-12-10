% Displays data for a single subject. This is a quick way to view the bad
% trials, bad electrodes and final PSDs for all protocols.

% badTrialRejectionFlag: 1: Don't reject badElectrodes, 2: reject badElectrodes for that protocol, 3: Reject badElectrodes of all protocols 

function displayPowerDataSingleSubject(subjectName,expDate,folderSourceString,badTrialNameStr,badElectrodeRejectionFlag,plotRawTFFlag,sortByBadTrialFlag)

if ~exist('folderSourceString','var');    folderSourceString=[];        end
if ~exist('badElectrodeList','var');      badTrialNameStr='_wo_v8';     end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=2;  end
if ~exist('plotRawTFFlag','var');         plotRawTFFlag=0;              end
if ~exist('sortByBadTrialFlag','var');    sortByBadTrialFlag=1;         end

if isempty(folderSourceString)
    folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG'; 
capType = 'actiCap64_UOL';

protocolNameList = [{'EO1'}     {'EC1'}     {'G1'}      {'M1'}          {'G2'}      {'EO2'}     {'EC2'}     {'M2'}];
colorNames       = [{[0.9 0 0]} {[0 0.9 0]} {[0 0 0.9]} {[0.9 0.9 0.9]} {[0 0 0.3]} {[0.3 0 0]} {[0 0.3 0]} {[0.3 0.3 0.3]}];
numProtocols = length(protocolNameList);

comparePSDConditions{1} = [1 6];
comparePSDConditions{2} = [2 7];
comparePSDConditions{3} = [3 5];
comparePSDConditions{4} = [4 8];
numPSDComparisons = length(comparePSDConditions);
comparePSDConditionStr = cell(1,numPSDComparisons);

for i=1:numPSDComparisons
   for j=1:length(comparePSDConditions{i})
       comparePSDConditionStr{i} = cat(2,comparePSDConditionStr{i},protocolNameList{comparePSDConditions{i}(j)});
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);
numGroups = length(electrodeGroupList);
electrodeGroupList{numGroups+1} = highPriorityElectrodeNums;
groupNameList{numGroups+1} = 'highPriority';
numGroups=numGroups+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hBadElectrodes = getPlotHandles(1,numProtocols,[0.05 0.875 0.6 0.1],0.01,0.01,1);
hBadElectrodes2 = getPlotHandles(1,4,[0.7 0.875 0.25 0.1],0.01,0.01,1);

hTF = getPlotHandles(numGroups,numProtocols,[0.05 0.05 0.6 0.8],0.01,0.01,1);
hPSD  = getPlotHandles(numGroups,numPSDComparisons,[0.7 0.05 0.25 0.8],0.01,0.01,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Ranges for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap jet;
freqRangeHz = [0 100];

if plotRawTFFlag
    cLims = [-3 3];
else
    cLims = [-1 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display bad electrodes for all protocols and also generate common bad electrodes
badTrialsList = cell(1,numProtocols);
badElecList = cell(1,numProtocols);
badElectrodes.badImpedanceElecs = [];
badElectrodes.noisyElecs = [];
badElectrodes.flatPSDElecs = [];
badElectrodes.declaredBadElecs = [];

for i=1:numProtocols
    protocolName=protocolNameList{i};
    badFileName = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData',['badTrials' badTrialNameStr '.mat']);
    if exist(badFileName,'file')
        x=getBadTrialInfo(badFileName);
        badTrialsList{i}=x.badTrials(:);
        badElecList{i} = x.badElecs;
        badElectrodes.badImpedanceElecs = cat(1,badElectrodes.badImpedanceElecs,x.badElecs.badImpedanceElecs);
        badElectrodes.noisyElecs = cat(1,badElectrodes.noisyElecs,x.badElecs.noisyElecs);
        badElectrodes.flatPSDElecs = cat(1,badElectrodes.flatPSDElecs,x.badElecs.flatPSDElecs);
        badElectrodes.declaredBadElecs = cat(1,badElectrodes.declaredBadElecs,x.badElecs.declaredBadElecs);
        displayBadElecs(hBadElectrodes(i),subjectName,expDate,protocolName,folderSourceString,gridType,capType,badTrialNameStr);
        title(hBadElectrodes(i),[protocolNameList{i} '(' num2str(length(getAllBadElecs(badElecList{i}))) ')'],'color',colorNames{i});
    else
        badTrialsList{i}=[];
        badElecList{i} = [];
    end
end

badElectrodes.badImpedanceElecs = unique(badElectrodes.badImpedanceElecs);
badElectrodes.noisyElecs = unique(badElectrodes.noisyElecs);
badElectrodes.flatPSDElecs = unique(badElectrodes.flatPSDElecs);
badElectrodes.declaredBadElecs = unique(badElectrodes.declaredBadElecs);
displayBadElecs(hBadElectrodes2(1),subjectName,expDate,protocolName,folderSourceString,gridType,capType,badTrialNameStr,badElectrodes,hBadElectrodes2(2));

% displayElectrodeGroups
showElectrodeGroups(hBadElectrodes2(3:4),capType,electrodeGroupList,groupNameList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Data
for g=1:numGroups

    disp(['Working on group: ' groupNameList{g}]);

    % Get Data
    psdVals = cell(1,numProtocols);
    freqVals = cell(1,numProtocols);
    meanPSDVals = cell(1,numProtocols);
    
    numGoodElectrodesList = zeros(1,numProtocols);
    numGoodTrials = zeros(1,numProtocols);

    for i=1:numProtocols
        
        if badElectrodeRejectionFlag==1
           electrodeList = electrodeGroupList{g};
        elseif badElectrodeRejectionFlag==2
            electrodeList = setdiff(electrodeGroupList{g},getAllBadElecs(badElecList{i}));
        elseif badElectrodeRejectionFlag==3
            electrodeList = setdiff(electrodeGroupList{g},getAllBadElecs(badElectrodes));
        end
        numGoodElectrodesList(i) = length(electrodeList);
        
        if ~isempty(electrodeList)
            protocolName = protocolNameList{i};
            [psdVals{i},freqVals{i}] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,electrodeList);
            goodPos = setdiff(1:size(psdVals{i},2),badTrialsList{i});
            meanPSDVals{i} = mean(psdVals{i}(:,goodPos),2);
            numGoodTrials(i) = length(goodPos);
        end
    end

    if isempty(meanPSDVals{1})
        disp('baseline condition does not exist. Plotting rawTF');
        plotRawTFFlag=1;
        cLims = [-3 3];
    end
    
    for i=1:numProtocols
        
        % Time-frequency plots
        if ~isempty(psdVals{i})
            numTrials = size(psdVals{i},2);

            if sortByBadTrialFlag
                xPos = [setdiff(1:numTrials,badTrialsList{i}) badTrialsList{i}'];
                badTrialIndicesToUse = (numTrials-length(badTrialsList{i})+1):numTrials;
            else
                xPos = 1:numTrials;
                badTrialIndicesToUse = badTrialsList{i};
            end

            if plotRawTFFlag
                pcolor(hTF(g,i),1:numTrials,freqVals{i},(psdVals{i}(:,xPos)));
            else
                bl = repmat(meanPSDVals{1},1,numTrials);
                pcolor(hTF(g,i),1:numTrials,freqVals{i},(psdVals{i}(:,xPos))-bl);
            end
            
            shading(hTF(g,i),'interp');
            clim(hTF(g,i),cLims); 
            ylim(hTF(g,i),freqRangeHz);
            
            hold(hTF(g,i),'on');
            plot(hTF(g,i),badTrialIndicesToUse,freqRangeHz(2)-1,'k.');
            text(1,freqRangeHz(2)-5,['Elec=' num2str(numGoodElectrodesList(i))],'parent',hTF(g,i));
        end

        if g==1
            title(hTF(g,i),['Good trials=' num2str(numGoodTrials(i))]);
        end

        if (i==1 && g<numGroups)
            set(hTF(g,i),'XTickLabel',[]); % only remove x label
        elseif (i>1 && g<numGroups)
            set(hTF(g,i),'XTickLabel',[],'YTickLabel',[]);
        elseif (i>1 && g==numGroups)
            set(hTF(g,i),'YTickLabel',[]); % only remove y label
        end
    end
    
    ylabel(hTF(g,1),groupNameList{g});
        
    % meanPSD Plots
    for i = 1:numPSDComparisons
        for j=1:length(comparePSDConditions{i})
            conditionNum = comparePSDConditions{i}(j);
            if ~isempty(meanPSDVals{conditionNum})
                if plotRawTFFlag
                    plot(hPSD(g,i),freqVals{conditionNum},meanPSDVals{conditionNum},'color',colorNames{conditionNum});
                else
                    plot(hPSD(g,i),freqVals{conditionNum},meanPSDVals{conditionNum}-meanPSDVals{1},'color',colorNames{conditionNum});
                end
            end
            hold(hPSD(g,i),'on');
        end
        
        if ~plotRawTFFlag
            plot(hPSD(g,i),freqVals{conditionNum},zeros(1,length(freqVals{conditionNum})),'k--');
        end
        
        xlim(hPSD(g,i),freqRangeHz);
        ylim(hPSD(g,i),cLims);
        if g<numGroups
            set(hPSD(g,i),'XTickLabel',[]);
        end
        
        if g==1
            title(hPSD(g,i),comparePSDConditionStr{i});
        end
    end
end

for i=1:numProtocols
    xlabel(hTF(numGroups,i),'TrialNum');
end
end

function x=getBadTrialInfo(badFileName)
x = load(badFileName);

if exist('getDeclaredBadElecs','file')
    x.badElecs.declaredBadElecs = getDeclaredBadElecs';
end
end
function [psd,freqVals] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,electrodeList)

timeRange = [0.25 1.25];
tapers = [1 1];
freqRange = [0 100];

folderExtract = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'extractedData');
folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');

if ~exist(folderExtract,'file')
    disp([folderExtract ' does not exist']);
    psd = []; freqVals=[];
else
    numElectrodes = length(electrodeList);
    
    t = load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    timeVals = t.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    goodTimePos = find(timeVals>=timeRange(1),1) + (1:round(Fs*diff(timeRange)));
    
    % Set up multitaper
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = freqRange;
    params.trialave = 0;
    
    for i=1:numElectrodes
        e = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeList(i)) '.mat']));
        [psdTMP(i,:,:),freqVals] = mtspectrumc(e.analogData(:,goodTimePos)',params); %#ok<AGROW>
    end
    psd = log10(squeeze(mean(psdTMP,1)));
end
end
function displayBadElecs(hBadElectrodes,subjectName,expDate,protocolName,folderSourceString,gridType,capType,badTrialNameStr,badElectrodes,hStats)

if ~exist('gridType','var');        gridType = 'EEG';                   end
if ~exist('capType','var');         capType = 'actiCap64_2019';         end
if ~exist('badTrialNameStr','var'); badTrialNameStr = '_v5';            end
if ~exist('badElectrodes','var');   badElectrodes = [];                 end
if ~exist('hStats','var');          hStats = [];                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');
badTrialsInfo = getBadTrialInfo(fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']));

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']); montageChanlocs = x.chanlocs;

if ~isfield(badTrialsInfo,'eegElectrodeLabels') % Check if the labels match with the save labels, if these labels have been saved
    disp('Electrode labels not specified in badTrials file. Taking from Montage...');
else
    if ~isequal(montageLabels(:),badTrialsInfo.eegElectrodeLabels(:))
        error('Montage labels do not match with channel labels in badTrials');
    else
%        disp('Montage labels match with channel labels in badTrials');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hBadElectrodes);
electrodeSize = 5;
numElectrodes = length(montageLabels);

for i=1:numElectrodes
%    label = num2str(i); %[num2str(i) '-' montageLabels{i}];
    montageChanlocs(i).labels = ' ';
end

% If you specify these, they override the existing ones
badImpedanceElectrodeColor = 'r';
noisyElectrodeColor = 'm';
flatPSDElectrodeColor = 'b';
declaredBadElectrodeColor = 'k';

if isempty(badElectrodes)
    badImpedanceElectrodes = badTrialsInfo.badElecs.badImpedanceElecs;
    noisyElectrodes = badTrialsInfo.badElecs.noisyElecs;
    flatPSDElectrodes = badTrialsInfo.badElecs.flatPSDElecs;
    declaredBadElectrodes = badTrialsInfo.badElecs.declaredBadElecs;
else
    badImpedanceElectrodes = badElectrodes.badImpedanceElecs;
    noisyElectrodes = badElectrodes.noisyElecs;
    flatPSDElectrodes = badElectrodes.flatPSDElecs;
    declaredBadElectrodes = badElectrodes.declaredBadElecs;
end

topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{badImpedanceElectrodes,'o',badImpedanceElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{noisyElectrodes,'o',noisyElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{flatPSDElectrodes,'o',flatPSDElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{declaredBadElectrodes,'o',declaredBadElectrodeColor,electrodeSize});
%topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{highPriorityElectrodeList,'o',highPriorityElectrodeColor,electrodeSize});
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

if ~isempty(hStats)
    axes(hStats)
    set(hStats,'visible','off');
    text(0.05,0.9,['badImpedance(' num2str(length(badImpedanceElectrodes)) ')'],'color',badImpedanceElectrodeColor);
    text(0.05,0.7,['Noisy(' num2str(length(noisyElectrodes)) ')'],'color',noisyElectrodeColor);
    text(0.05,0.5,['FlatPSD(' num2str(length(flatPSDElectrodes)) ')'],'color',flatPSDElectrodeColor);
    text(0.05,0.3,['declaredBad(' num2str(length(declaredBadElectrodes)) ')'],'color',declaredBadElectrodeColor);
    text(0.05,0.1,['All(' num2str(length(getAllBadElecs(badElectrodes))) ')'],'color',declaredBadElectrodeColor);
end
end
function allBadElecs = getAllBadElecs(badElectrodes)
if ~isempty(badElectrodes)
    allBadElecs = unique([badElectrodes.badImpedanceElecs; badElectrodes.noisyElecs; badElectrodes.flatPSDElecs; badElectrodes.declaredBadElecs]);
else
    allBadElecs = [];
end
end
function showElectrodeGroups(hPlots,capType,electrodeGroupList,groupNameList)

if ~exist('capType','var');         capType = 'actiCap64_2019';         end

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
    topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize});
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

axes(hPlots(2))
set(hPlots(2),'visible','off');
for i=1:numElectrodeGroups
    text(0.05,0.9-0.15*(i-1),groupNameList{i},'color',electrodeGroupColorList(i,:),'unit','normalized');
end

end