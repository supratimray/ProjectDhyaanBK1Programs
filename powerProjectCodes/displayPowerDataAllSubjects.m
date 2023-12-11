% analysisChoice - 'st', 'bl' or 'combined'

% refChoice - 'none' (show raw PSDs and power) or a protocolName. For
% example, if refChoice is 'G1' then we use G1 baseline as reference.

% badTrialRejectionFlag: 
% 1: Reject badElectrodes of protocolName
% 2. Reject common badElectrodes of all protocols
% 3: Reject badElectrodes of G1

function displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots)

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('refChoice','var');             refChoice='none';             end

if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=1;  end

if ~exist('stRange','var');               stRange = [0.25 1.25];        end

if ~exist('freqRangeList','var')    
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [22 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [5 50];
end
cutoffNumElectrodes = cutoffList(1);
cutoffNumTrials = cutoffList(2);

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
protocolPos = find(strcmp(protocolNameList,protocolName));

if ~strcmp(refChoice,'none')
    protocolPosRef = find(strcmp(protocolNameList,refChoice));
else
    protocolPosRef = [];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        if ~isempty(protocolPosRef)
            tmpPowerRef = getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
        end

        if isempty(protocolPosRef) % No need to worry about Ref
            if isempty(tmpPower)
                disp(['Not enough trials for subject: ' subjectName]);
            else
                powerDataTMP = cat(3,powerDataTMP,tmpPower);
            end
        else
            if isempty(tmpPowerRef) || isempty(tmpPower) % If either one is empty
                disp(['Not enough trials for protocol or ref condition of subject: ' subjectName]);
            else
                powerDataTMP = cat(3,powerDataTMP,tmpPower);
                powerDataRefTMP = cat(3,powerDataRefTMP,tmpPowerRef);
            end
        end
    end
    powerData{i} = powerDataTMP;
    powerDataRef{i} = powerDataRefTMP;
end

%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
for i = 1:numFreqRanges
    freqPosList{i} = intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElectrodes = size(powerData{1},1);
percentData = zeros(2,numElectrodes);
comparisonData = zeros(numFreqRanges,2,numElectrodes);

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
    axes(hTopo1(i)); %#ok<*LAXES>

    % Modification in the topoplot code which allows us to not interpolate across electrodes.
    %topoplot_murty(numBadSubjects/numSubjects,montageChanlocs,'electrodes','off','style','blank','drawaxis','off','emarkercolors',numBadSubjects/numSubjects); colorbar;
    percentData(i,:) = 100*(numBadSubjects/numSubjects);
    topoplot(percentData(i,:),montageChanlocs,'maplimits',[0 100],'electrodes','on'); colorbar;
    title(titleStr{i},'color',displaySettings.colorNames(i,:));
    if i==1
        ylabel('Bad subjects (%)');
    end

    %%%%%%%%%%%%%%%%%%%%%%% Show topoplots of power %%%%%%%%%%%%%%%%%%%%%%%
    for j=1:numFreqRanges
        axes(hTopo2(j,i));
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
        topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo); colorbar;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plot the difference of topoplots %%%%%%%%%%%%%%%%%%%
axes(hTopo1(3));
topoplot(-diff(percentData),montageChanlocs,'maplimits',[-25 25],'electrodes','on'); colorbar;

for i=1:numFreqRanges
    axes(hTopo2(i,3));
    data = -diff(squeeze(comparisonData(i,:,:)));
    if isempty(protocolPosRef)
        data = 10*data;
    end
    topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo); colorbar;
end

%%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numGroups
    meanPSDData = cell(1,2);
    meanPSDDataRef = cell(1,2);
    logPSDData = cell(1,2);
    
    for j=1:2
        pData = powerData{j}(electrodeGroupList{i},:,:);
        numGoodElecs = length(electrodeGroupList{i}) - sum(isnan(squeeze(pData(:,1,:))),1);
        badSubjectPos = find(numGoodElecs<=cutoffNumElectrodes);
        
        if ~isempty(protocolPosRef)
            pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
            numGoodElecsRef = length(electrodeGroupList{i}) - sum(isnan(squeeze(pDataRef(:,1,:))),1);
            badSubjectPosRef = find(numGoodElecsRef<=cutoffNumElectrodes);
            badSubjectPos = unique(cat(2,badSubjectPos,badSubjectPosRef));
        end

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

        text(30,yLimsPSD(2)-0.5*j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(i));
    end
    displayAndcompareData(hPSD(i),logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,1);
    title(groupNameList{i});
    xlim(hPSD(i),freqLims);

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

        % display violin plots for power
        displaySettings.plotAxes = hPower(j,i);
        if i==numGroups && j==1
            displaySettings.showYTicks=1;
            displaySettings.showXTicks=1;
        else
            displaySettings.showYTicks=0;
            displaySettings.showXTicks=0;
        end
        displayViolinPlot(tmpLogPower,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,0,displaySettings);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize});
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

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
function displayAndcompareData(hPlot,data,xs,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,nonMatchedFlag)

if ~exist('displaySignificanceFlag','var'); displaySignificanceFlag=0;  end
if ~exist('useMedianFlag','var');           useMedianFlag=1;            end
if ~exist('nonMatchedFlag','var');          nonMatchedFlag=1;           end

if useMedianFlag
    getLoc = @(g)(squeeze(median(g,1)));
else
    getLoc = @(g)(squeeze(mean(g,1)));
end

numGroups = length(data);

axes(hPlot);
for i=1:numGroups
    clear bootStat mData sData
    mData = getLoc(data{i}); 
    if useMedianFlag
        bootStat = bootstrp(1000,getLoc,data{i});
        sData = std(bootStat);
    else
        sData = std(data{i},[],1)/sqrt(size(data{i},1));
    end
    
    patch([xs';flipud(xs')],[mData'-sData';flipud(mData'+sData')],displaySettings.colorNames(i,:),'linestyle','none','FaceAlpha',0.4);
    hold on;
    plot(xs,mData,'color',displaySettings.colorNames(i,:),'linewidth',1);
end

set(gca,'fontsize',displaySettings.fontSizeLarge);
set(gca,'TickDir','out','TickLength',displaySettings.tickLengthMedium);

if exist('yLims','var') && ~isempty(yLims)
    ylim(yLims);
else
    yLims = ylim;
end

if displaySignificanceFlag % Do significance Testing
    
    allData = [];
    allIDs = [];
    for j=1:numGroups
        allData = cat(1,allData,data{j});
        allIDs = cat(1,allIDs,j+zeros(size(data{j},1),1));
    end
       
   for i=1:length(xs)
       if useMedianFlag
           p=kruskalwallis(allData(:,i),allIDs,'off');
       else
           if nonMatchedFlag
               [~,p]=ttest2(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
           else
               [~,p]=ttest(data{1}(:,i),data{2}(:,i)); % only tests 2 groups
           end
       end
       % Get patch coordinates
       yVals = yLims(1)+[0 0 diff(yLims)/20 diff(yLims)/20];
       
       clear xMidPos xBegPos xEndPos
       xMidPos = xs(i);
       if i==1
           xBegPos = xMidPos;
       else
           xBegPos = xMidPos-(xs(i)-xs(i-1))/2; 
       end
       if i==length(xs)
           xEndPos = xMidPos; 
       else
           xEndPos = xMidPos+(xs(i+1)-xs(i))/2; 
       end
       clear xVals; xVals = [xBegPos xEndPos xEndPos xBegPos]';
       
       if (p<0.05)
           patch(xVals,yVals,'c','linestyle','none');
       end
       if (p<0.01)
           patch(xVals,yVals,'k','linestyle','none');
       end
   end
end
end