% analysisChoice - 'st', 'bl' or 'combined'
% badTrialRejectionFlag: 
% 1: Don't reject badElectrodes
% 2: Reject badElectrodes of protocolName
% 3. Reject common badElectrodes of all protocols
% 4: Reject badElectrodes of G1

function displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,useMedianFlag)

if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('badEyeCondition','var');       badEyeCondition='wo';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=2;  end
if ~exist('stRange','var');               stRange = [0.25 1.25];        end
if ~exist('freqRangeList','var')       
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 35]; % Slow gamma (SG)
    freqRangeList{3} = [40 70]; % Fast gamma (FG)
    freqRangeList{4} = [80 150]; % high-gamma (HG)
end
if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end

numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10; 
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [1 0 0];
displaySettings.colorNames(2,:) = [0 1 0];
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';
cLimsTopo = [-1.5 1.5];
yLimsPSD = [-2.5 2.5];

%%%%%%%%%%%%%%%%%%%%%%%% Get electrode groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';
saveFolderName = 'savedData';

[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);
numGroups = length(electrodeGroupList);
electrodeGroupList{numGroups+1} = highPriorityElectrodeNums;
groupNameList{numGroups+1} = 'highPriority';
numGroups=numGroups+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPSD  = getPlotHandles(numGroups,1,[0.05 0.05 0.25 0.9],0.02,0.02,1);
hPower = getPlotHandles(numGroups,numFreqRanges,[0.325 0.05 0.35 0.9],0.02,0.02,1);
hTopo = getPlotHandles(numFreqRanges+2,2,[0.7 0.05 0.25 0.9],0.02,0.02,1); % First two plots show the layout and fraction of good subjects per electrode

montageChanlocs = showElectrodeGroups(hTopo(1,:),capType,electrodeGroupList,groupNameList);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocol Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
protocolPos = find(strcmp(protocolNameList,protocolName));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutoffNumTrials = 50;

powerData = cell(1,2);
for i=1:2
    powerDataTMP=[];

    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));
        numTrials = tmpData.numTrials(protocolPos);
        freqVals = tmpData.freqVals;
        badElectrodes = getBadElectrodes(tmpData.badElectrodes,badElectrodeRejectionFlag,protocolPos);

        if numTrials < cutoffNumTrials
            disp(['Not enough trials for subject: ' subjectName]);
        else
            if strcmpi(analysisChoice,'st')
                tmpPower = tmpData.psdValsST{protocolPos};
            elseif strcmpi(analysisChoice,'bl')
                tmpPower = tmpData.psdValsBL{protocolPos};
            else
                tmpPower = (tmpData.psdValsST{protocolPos}+tmpData.psdValsBL{protocolPos})/2; % average
            end
            tmpPower(badElectrodes,:) = NaN;
            powerDataTMP = cat(3,powerDataTMP,tmpPower);
        end
    end
    powerData{i} = powerDataTMP;
end

%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
for i = 1:numFreqRanges
    freqPosList{i} = intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2
    x=powerData{i};
    numSubjects = size(x,3);
    numElectrodes = size(x,1);

    %%%%%%%%%%%% Show percent of bad subjects per electrode %%%%%%%%%%%%%%%
    numBadSubjects = zeros(1,numElectrodes);
    for j=1:numElectrodes
        numBadSubjects(j) = sum(isnan(squeeze(x(j,1,:))));
    end
    axes(hTopo(2,i)); %#ok<*LAXES>
    % Modification in the topoplot code which allows us to not interpolate across electrodes.
    %topoplot_murty(numBadSubjects/numSubjects,montageChanlocs,'electrodes','off','style','blank','drawaxis','off','emarkercolors',numBadSubjects/numSubjects); colorbar;
    topoplot(100*(numBadSubjects/numSubjects),montageChanlocs,'maplimits',[0 100],'electrodes','on'); colorbar;
    title(titleStr{i},'color',displaySettings.colorNames(i,:));
    if i==1
        ylabel('Bad subjects (%)');
    end

    %%%%%%%%%%%%%%%%%%%%%%% Show topoplots of power %%%%%%%%%%%%%%%%%%%%%%%
    for j=1:numFreqRanges
        axes(hTopo(2+j,i));
        if useMedianFlag
            data = squeeze(median(log10(sum(powerData{i}(:,freqPosList{j},:),2)),3,'omitnan'));
        else
            data = squeeze(mean(log10(sum(powerData{i}(:,freqPosList{j},:),2)),3,'omitnan'));
        end
        topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo); colorbar;
    end
end

%%%%%%%%%%%%%%%%%%%%%% Plots PSDs and power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutoffNumElectrodes = 3;

for i=1:numGroups
    meanPSDData = cell(1,2);
    logPSDData = cell(1,2);
    
    for j=1:2
        pData = powerData{j}(electrodeGroupList{i},:,:);
        numGoodElecs = length(electrodeGroupList{i}) - sum(isnan(squeeze(pData(:,1,:))),1);
        badSubjectPos = find(numGoodElecs<=cutoffNumElectrodes);
        
        if ~isempty(badSubjectPos)
            disp([groupNameList{i} ', ' titleStr{j} ', '  'Not enough good electrodes for ' num2str(length(badSubjectPos)) ' subjects.']);
            pData(:,:,badSubjectPos)=[];
        end
        meanPSDData{j} = squeeze(mean(pData,1,'omitnan'))';
        logPSDData{j} = log10(meanPSDData{j});

        text(freqVals(end)-50,yLimsPSD(2)-j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(i));
    end
    displayAndcompareData(hPSD(i),logPSDData,freqVals,displaySettings,yLimsPSD,1,useMedianFlag,1);
    ylabel(groupNameList{i});

    % Violin plots for power 
    for j=1:numFreqRanges
        tmpLogPower = cell(1,2);
        for k=1:2
            tmpLogPower{k} = log10(squeeze(sum(meanPSDData{k}(:,freqPosList{j}),2)));
        end
        axes(hPower(i,j));
        displayViolinPlot(tmpLogPower,[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,0);
        if i==1
            title(hPower(i,j),[num2str(freqRangeList{j}(1)) '-' num2str(freqRangeList{j}(2)) ' Hz'],'color',freqRangeColors(j,:));
        end

        % Add lines in PSD plots
        for k=1:2
            line([freqRangeList{j}(k) freqRangeList{j}(k)],yLimsPSD,'color',freqRangeColors(j,:),'parent',hPSD(i));
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if badElectrodeRejectionFlag==1 % No bad electrodes
    badElectrodes = [];
elseif badElectrodeRejectionFlag==2 % Bad electrodes for the protocol
    badElectrodes = badElectrodeList{protocolPos};
elseif badElectrodeRejectionFlag==3 % common bad electrodes for all protocols
    badElectrodes=[];
    for i=1:length(badElectrodeList)
        badElectrodes=cat(1,badElectrodes,badElectrodeList{i});
    end
    badElectrodes = unique(badElectrodes);
elseif badElectrodeRejectionFlag==4 % Bad electrodes of G1
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
           patch(xVals,yVals,'k','linestyle','none');
       end
       if (p<0.01)
           patch(xVals,yVals,'g','linestyle','none');
       end
   end
end
end