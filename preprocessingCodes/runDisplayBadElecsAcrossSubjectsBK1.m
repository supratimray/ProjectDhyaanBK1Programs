% runDisplayBadElecsAcrossSubjectsBK1
% This script displays bad electrodes across subjects
% Displays for the paired subject if the 'dispForPairedSubject' flag is 'ON'; default: all the subjects
% Can sort the subjects according to date of data collection if 'sortByDate' flag is 'ON';
% default: displays for the Meditators followed by Controls
% Also, displays the badSubjects and bad electrodes according to the threshold value
% Added functionality to display the declaredBadElecs based on the threshold value (set removeDeclaredBadElecs flag=0)
% For paired subject, if we declare bad electrodes the difference is not significant

close all
clear
fh = figure(1);
fh.WindowState = 'maximized';
% dispForPairedSubject = 0; removeDeclaredBadElecs = 0; removeDeclaredBadSubjects = 0; % Use this option to get all 76 subjects and all electrodes
% dispForPairedSubject = 0; removeDeclaredBadElecs = 1; removeDeclaredBadSubjects = 0; % Use this option to get all 76 subjects but remove the electrodes that are declared to be bad
dispForPairedSubject = 0; removeDeclaredBadElecs = 1; removeDeclaredBadSubjects = 1; % Use this option to get 74 subjects and remove the electrodes that are declared to be bad

sortByDate = 0;
badElecThresoldBadSub = 0.50;
badElecThresoldAcrossSubject = 0.35;

% fixed variables
badTrialNameStr = '_wo_v8';
numElecs = 64;
allEEGElecArray = 1:numElecs;
protocolName = 'G1';
gridType = 'EEG';
capType = 'actiCap64_UOL';
folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
[~,~,~,electrodeGroupList0,groupNameList0,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType);
numElecGroups =  length(electrodeGroupList0);
[subjectNameList,expDateList] = getDemographicDetails('BK1');

% get the subject List
if dispForPairedSubject % paired
    pairedSubjectNameList = getPairedSubjectsBK1;
    goodSubjectList = pairedSubjectNameList(:);
    controlSubjectList = pairedSubjectNameList(:,2);
    numSubjects = length(goodSubjectList);
    numMeditators = size(pairedSubjectNameList,1);
    axisTitle = sgtitle(['Paired Subjects: n = ' num2str(numSubjects)]);
    set(axisTitle, 'FontSize', 20);
else % all the subjects
    if removeDeclaredBadSubjects
        [allSubjectList, meditatorList, controlList] = getGoodSubjectsBK1;
    else
        fileName = 'BK1AllSubjectList.mat';
        load(fileName,'allSubjectList','controlList','meditatorList');
    end
    controlSubjectList = controlList;
    goodSubjectList = allSubjectList;
    numSubjects = length(goodSubjectList);
    numMeditators = length(meditatorList);
    axisTitle = sgtitle(['All Subjects: n = ' num2str(numSubjects)]);
    set(axisTitle, 'FontSize', 20);
end
numControls = numSubjects-numMeditators;
% BadElecs x Subject matrix
%-----------------------------------------------------
absExpDate        = zeros(1,numSubjects);
badSubjectStatus  = zeros(1,numSubjects);
badElecPercentage = zeros(1,numSubjects);
allBadElecsMatrix = zeros(numElecs,numSubjects);
for s=1:numSubjects
    subjectName = goodSubjectList{s};
    disp(['Extracting Subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,subjectNameList)};
    folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');
    badTrialsInfo = load(fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']));
    % bad elecs:
    noisyElecs        = badTrialsInfo.badElecs.noisyElecs;
    flatPSDElecs      = badTrialsInfo.badElecs.flatPSDElecs;
    badImpedanceElecs = badTrialsInfo.badElecs.badImpedanceElecs;
    if removeDeclaredBadElecs
        declaredBadElectrodes = getDeclaredBadElecs;
        badElecrodesIndex     = unique([badImpedanceElecs;noisyElecs;flatPSDElecs;declaredBadElectrodes']);
    else
        badElecrodesIndex = unique([badImpedanceElecs;noisyElecs;flatPSDElecs]);
    end

    % assign the values to the matrix
    allBadElecsMatrix(badElecrodesIndex,s) = 1;
    % calculate and save bad elec percentage
    % and assign the badSubject Status
    badElecPercentage(1,s) = length(badElecrodesIndex)/numElecs;
    if badElecPercentage(1,s) > badElecThresoldBadSub
        badSubjectStatus(1,s) = 1;
    end
    % save the absolute date of recording
    x = expDate;
    absExpDate(1,s) = datenum([x(1:2) '/' x(3:4) '/' x(5:6)],'dd/mm/yy');
end

% get the bad elecs percentage for individual subjects per group
badElecPercentageIndMeditators = badElecPercentage(:,1:numMeditators);
badElecPercentageIndConrols    = badElecPercentage(:,numMeditators+1:end);

% gets the bad elecs percentage across subjects
% badElecPercentageAcrossSubjects       = round(sum(allBadElecsMatrix,2)/numSubjects,2);
% binarybadElecPercentageAcrossSubjects = (badElecPercentageAcrossSubjects > badElecThresoldAcrossSubject);

% get the bad elecs percentage across groups
badElecPercentageAcrossMeditators       = round(sum(allBadElecsMatrix(:,1:numMeditators),2)/numMeditators,2);
badElecPercentageAcrossControls         = round(sum(allBadElecsMatrix(:,numMeditators+1:end),2)/numControls,2);
binarybadElecPercentageAcrossSubjects   = (badElecPercentageAcrossMeditators > badElecThresoldAcrossSubject) | (badElecPercentageAcrossControls > badElecThresoldAcrossSubject);
binaryDeclaredBadElecAcrossSubjects     = ismember(allEEGElecArray, declaredBadElectrodes)';

badElecsCurrentThreshold                = allEEGElecArray(binarybadElecPercentageAcrossSubjects);
disp(['Bad elctrodes according to the current thresold(' num2str(badElecThresoldAcrossSubject) ') are ' num2str(badElecsCurrentThreshold)]);

% sorting the subjects accoring to exp dates
[sortedExpDate,orginalSortedIndexExpDate] = sort(absExpDate);
sortedbadElecPercent   = badElecPercentage(orginalSortedIndexExpDate);
sortedBadSubjectStatus = badSubjectStatus(orginalSortedIndexExpDate);
sortedSubjectList      = goodSubjectList(orginalSortedIndexExpDate);

numExpDate=numSubjects;
subjectGroupMatrix = zeros(1,numSubjects);
allBadElecsMatrixSortedExpDate = zeros(numElecs,numSubjects);
for e=1:numExpDate
    if sortByDate
        allBadElecsMatrixSortedExpDate(:,e) = allBadElecsMatrix(:,orginalSortedIndexExpDate(e));
        subjectName = sortedSubjectList(e);
    else
        allBadElecsMatrixSortedExpDate(:,e)=allBadElecsMatrix(:,e);
        subjectName = goodSubjectList(e);
    end
    % assign the control subjects 1; Meditators by default is 0
    if any(strcmp(controlSubjectList,subjectName))
        subjectGroupMatrix(1,e) = 1;
    end
end

% sorting the electrodes according to the different scalp areas
elecVal = 1;
startPos = 1;
elecGroupVec = zeros(numElecs,1);
sortedElecMatrix = zeros(numElecs,numSubjects);
% badDeclaredElecStatus = zeros(numElecs,numSubjects);

for g=1:numElecGroups % Electrode Group
    electrodeList = electrodeGroupList0{g};
    endPos = startPos+length(electrodeList)-1;

    sortedElecMatrix(startPos:endPos,:)         = allBadElecsMatrixSortedExpDate(electrodeList,:);
    badElecStatus(startPos:endPos,:)            = binarybadElecPercentageAcrossSubjects(electrodeList,:);
    badDeclaredElecStatus(startPos:endPos,:)    = binaryDeclaredBadElecAcrossSubjects(electrodeList,:);

    elecGroupVec(startPos:endPos,:) = elecVal;
    startPos = startPos+length(electrodeList);
    elecVal = elecVal+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot elcs x Sub Matrix
h1 = getPlotHandles(1,1,[0.15 0.1 0.8 0.7]);
imagesc(1:numSubjects,numElecs:-1:1,flipud(sortedElecMatrix),'parent',h1);
set(h1,'Xtick',1:1:numSubjects,'TickDir','out','TickLength',[0.005, 0.001],'XTickLabelRotation',0);
set(h1,'Ytick',1:1:numElecs,'YTickLabel',[],'TickDir','out','TickLength',[0.005, 0.001]);

colormap(h1,[gray(64)]);
set(gca,'YDir','normal');
if sortByDate
    xlabel('Subject Index (Sorted According to experiment date)','fontsize',15,'fontweight','bold');
else
    xlabel('Subject Index (Sorted According to the subject category)','fontsize',15,'fontweight','bold');
end

% shows the badSubjects
yPos = 71;
xOffSet = 0.265;
if sortByDate
    badSubInd = find(sortedBadSubjectStatus==1);
    for b=1:length(badSubInd)
        text(badSubInd(b)-xOffSet,yPos,['x (' num2str(sortedbadElecPercent(badSubInd(b))) ')'],'Color','Red','FontSize',12,'FontWeight','bold','Parent',h1);
        text(badSubInd(b)-xOffSet,yPos+2,sortedSubjectList{badSubInd(b)},'Color','Red','FontSize',12,'FontWeight','bold','Parent',h1);
    end
else
    badSubInd = find(badSubjectStatus==1);
    for b=1:length(badSubInd)
        text(badSubInd(b)-xOffSet,yPos,['x (' num2str(badElecPercentage(badSubInd(b))) ')'],'Color','Red','FontSize',12,'FontWeight','bold','Parent',h1);
        text(badSubInd(b)-xOffSet,yPos+2,goodSubjectList{badSubInd(b)},'Color','Red','FontSize',12,'FontWeight','bold','Parent',h1);
    end
end


% shows the badElectrodes
xPos=numSubjects+1;
yOffSet = 0.4;
badElecInd = find(badElecStatus==1);
declaredBadElecInd = find(badDeclaredElecStatus==1);
newBadElecInd = setdiff(badElecInd,declaredBadElecInd);
for e=1:length(badElecInd)
    if ismember(badElecInd(e),newBadElecInd)
        text(xPos,badElecInd(e)+yOffSet,'x','Color','Green','FontSize',12,'FontWeight','bold','Parent',h1);
    else
        text(xPos,badElecInd(e)+yOffSet,'x','Color','Red','FontSize',12,'FontWeight','bold','Parent',h1);
    end 
end

% shows the elecGroups
h0 = getPlotHandles(1,1,[0.11 0.1 0.03 0.7]);
subplot(h0);
imagesc(1:numSubjects,length(elecGroupVec):-1:1,flipud(elecGroupVec),'parent',h0);
colormap(h0,[turbo(64)]);

yPos = 65;
xPos = -540;
delta = 15;
for g=1:numElecGroups
    groupName = groupNameList0{g};
    text(xPos,yPos,groupName,'FontSize',14,'Rotation',45,'parent',h0);
    if g==4
        delta = 10;
    end
    yPos = yPos-delta;
end
set(h0,'xTick',[],'Xticklabel',[],'YTick',[],'Yticklabel',[]);

% shows subject groups
h2 = getPlotHandles(1,1,[0.15 0.05+0.75 0.8 0.05]);
imagesc(1:numSubjects,length(allBadElecsMatrix):-1:1,flipud(subjectGroupMatrix),'parent',h2);
set(h2,'Xtick',1:1:numSubjects,'XTickLabel',[],'TickDir','out','TickLength',[0.005, 0.001]);

% Make a custom colorMap
colorMapSubject = zeros(64,3);
colorMapSubject(1:32,:)  = repmat([1 0 0],32,1);
colorMapSubject(33:64,:) = repmat([0 1 0],32,1);

colormap(h2,colorMapSubject);
set(gca,'YDir','normal','YTick',[],'Yticklabel',[]);

% denotes the subject Groups
text(5,300,'Meditators','Color','red','Rotation',0,'FontSize',14,'fontweight','bold','parent',h2);
text(5,410,'Controls','Color','green','Rotation',0,'FontSize',14,'fontweight','bold','parent',h2);

%%%%%%%%%%%%%%%% Display the bad electrodes across groups %%%%%%%%%%%%%%%%%
dataArray{1,1} = badElecPercentageIndMeditators * numElecs;
dataArray{1,2} = badElecPercentageIndConrols * numElecs;

colorArray = {'r','g'};
showData = 1;
plotQuartiles = 1;
showSignificance = 1;
pairedDataFlag = dispForPairedSubject;

fh2 = figure(2);
fh2.WindowState = 'maximized';
figHandles=displayViolinPlot(dataArray,colorArray,showData,plotQuartiles,showSignificance,pairedDataFlag);
ylabel('Num Bad Elecs');
xticks(1:length(dataArray));
xticklabels({'Meditators', 'Controls'});
set(gca, 'FontSize', 12,'FontWeight','Bold');
set(gca, 'TickDir', 'out');