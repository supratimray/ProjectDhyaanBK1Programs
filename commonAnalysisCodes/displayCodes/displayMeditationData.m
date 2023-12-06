function displayMeditationData(subjectName,expDate,folderSourceString,badTrialNameStr,badElectrodeRejectionFlag,rawTfTrialWise,refScheme,trialAvgFlag,diffTf)

if ~exist('folderSourceString','var');              folderSourceString = [];                          end
if ~exist('badTrialNameStr','var');                 badTrialNameStr = '_v5';                          end
if ~exist('badElectrodeRejectionFlag','var');       badElectrodeRejectionFlag = 2;                    end
if ~exist('plotRawTFFlag','var');                   rawTfTrialWise = 0;                               end
if ~exist('diffTf','var');                          diffTf = 0;                                       end
if ~exist('refScheme','var');                       refScheme = 1;                                    end
if ~exist('trialAvgFlag','var');                    trialAvgFlag = 0;                                 end

if isempty(folderSourceString)
    folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\MeditationProjects\MeditationProject2';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fixed variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';

if refScheme == 2
    load(['bipChInfo' capType '.mat']);
    allElectrodeList = 1:length(bipolarLocs);
else
    allElectrodeList = 1:64; % All EEG electrodes
end

protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
colorNames       = [{[0.9 0 0]} {[0 0.9 0]} {[0 0 0.9]} {[0.7 0.7 0.7]} {[0 0 0.3]} {[0.3 0 0]} {[0 0.3 0]} {[0.3 0.3 0.3]}];

numProtocols = length(protocolNameList);

% PSD comparisons between segments

comparePSDConditions{1} = [1 6];
comparePSDConditions{2} = [2 7];
comparePSDConditions{3} = [3 5];
comparePSDConditions{4} = [4 8];
numPSDComparisons = length(comparePSDConditions);
comparePSDConditionStr = cell(1,numPSDComparisons);

for i=1:numPSDComparisons
    for s=1:length(comparePSDConditions{i})
        comparePSDConditionStr{i} = cat(2,comparePSDConditionStr{i},protocolNameList{comparePSDConditions{i}(s)});
    end
end

freqList{1} = [8 12];  freqListNames{1} = 'Alpha';
freqList{2} = [20 34]; freqListNames{2} = 'SG';
freqList{3} = [36 66]; freqListNames{3} = 'FG';
numFreqRanges = length(freqList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType,refScheme);
numGroups = length(electrodeGroupList);
electrodeGroupList{numGroups+1} = highPriorityElectrodeNums;
groupNameList{numGroups+1} = 'highPriority';
numGroups = numGroups+1;
subSegments = 4; % no of subSements (M1 and M2 are divided into M1a, M1b and M1c)

%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hBadElectrodes  = getPlotHandles(1,numProtocols+subSegments,[0.05 0.875 0.6 0.1],0.01,0.01,1);
hBadElectrodes2 = getPlotHandles(1,4,[0.7 0.875 0.25 0.1],0.01,0.01,1);

hTF   = getPlotHandles(numGroups,numProtocols+subSegments,[0.05 0.35 0.6 0.5],0.01,0.01,1);
hPSD  = getPlotHandles(numGroups,numPSDComparisons,[0.7 0.35 0.25 0.5],0.01,0.01,1);
hTopo = getPlotHandles(numFreqRanges,numProtocols+subSegments,[0.05 0.05 0.6 0.25],0.01,0.01,1);
hPowerVsTime = getPlotHandles(numFreqRanges,1,[0.7 0.05 0.25 0.25],0.01,0.01,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Ranges for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colormap jet;
freqRangeHz = [0 100];

if rawTfTrialWise
    cLimsPlot = [-3 3];
else
    cLimsPlot = [-5 5];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display bad electrodes for all protocols and also generate common bad electrodes
badTrialsList = cell(1,numProtocols);
badElecList   = cell(1,numProtocols);
badElectrodes.badImpedanceElecs = [];
badElectrodes.noisyElecs        = [];
badElectrodes.flatPSDElecs      = [];

for i=1:numProtocols
    protocolName = protocolNameList{i};
    badFileName  = fullfile(folderSourceString,'data',subjectName,gridType, ...
        expDate,protocolName,'segmentedData',['badTrials' badTrialNameStr '.mat']);
    % display bad electrodes
    if exist(badFileName,'file')
        x = load(badFileName);
        badTrialsList{i} = x.badTrials;
        badElecList{i}   = x.badElecs;
        badElectrodes.noisyElecs        = cat(1,badElectrodes.noisyElecs,x.badElecs.noisyElecs);
        badElectrodes.flatPSDElecs      = cat(1,badElectrodes.flatPSDElecs,x.badElecs.flatPSDElecs);
        badElectrodes.badImpedanceElecs = cat(1,badElectrodes.badImpedanceElecs,x.badElecs.badImpedanceElecs);

        if i>4
            j=i+2; % skip two segments M1.b and M1.c to G2
            set(hBadElectrodes(5), 'visible','off');
            set(hBadElectrodes(6), 'visible','off');
            set(hBadElectrodes(11),'visible','off');
            set(hBadElectrodes(12),'visible','off');
        else
            j=i;
        end
        displayBadElecs(hBadElectrodes(j),subjectName,expDate,protocolName, ...
            folderSourceString,gridType,capType,badTrialNameStr);
        title(hBadElectrodes(j),protocolNameList{i},'color',colorNames{i});
    else
        badTrialsList{i} = [];
        badElecList{i}   = [];
    end
end
badElectrodes.badImpedanceElecs = unique(badElectrodes.badImpedanceElecs);
badElectrodes.noisyElecs        = unique(badElectrodes.noisyElecs);
badElectrodes.flatPSDElecs      = unique(badElectrodes.flatPSDElecs);
badElectrodes.declaredBadElecs  =  getDeclaredBadElecs(0.35);
displayBadElecs(hBadElectrodes2(1),subjectName,expDate,protocolName, ...
    folderSourceString,gridType,capType,badTrialNameStr,badElectrodes,hBadElectrodes2(2));

% displayElectrodeGroups
showElectrodeGroups(hBadElectrodes2(3:4),capType,electrodeGroupList,groupNameList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting the Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for g=1:numGroups % for different electrode groups getting the power data

    if g==4
        flag;
    end

    disp(['Working on group: ' groupNameList{g}]);
    psdVals     = cell(1,numProtocols);
    freqVals    = cell(1,numProtocols);
    meanPSDVals = cell(1,numProtocols);

    meanAlphaCurrSeg = zeros(8,5);
    meanSGCurrSeg    = zeros(8,5);
    meanFGCurrSeg    = zeros(8,5);

    numGoodElectrodesList = zeros(1,numProtocols);

    for i=1:numProtocols
        % reject bad electrodes
        if badElectrodeRejectionFlag==1
            electrodeList = electrodeGroupList{g};
        elseif badElectrodeRejectionFlag==2
            electrodeList = setdiff(electrodeGroupList{g},getAllBadElecs(badElecList{i}));
        elseif badElectrodeRejectionFlag==3
            electrodeList = setdiff(electrodeGroupList{g},getAllBadElecs(badElectrodes));
        elseif badElectrodeRejectionFlag==4
            % reject declared bad electrodes as well as protocol specific bad electrodes
            declaredBadElectrodes = getDeclaredBadElecs(0.35);
            electrodeList = setdiff(electrodeGroupList{g},union(declaredBadElectrodes,getAllBadElecs(badElecList{i})));
        end

        numGoodElectrodesList(i) = length(electrodeList);

        if ~isempty(electrodeList)
            protocolName = protocolNameList{i};
            badTrials = badTrialsList{i};
            disp(protocolName);

            [psdVals{i},freqVals{i},psdAcrossElc,tfPower{i},timeValsTF,~] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,electrodeList,refScheme,capType,trialAvgFlag,diffTf,badTrials);

            meanPSDVals{i} = mean(psdVals{i}(:,setdiff(1:size(psdVals{i},2),badTrialsList{i})),2);      % psd across good trials

            if g==6
                % For the highPriority electrode set, segmenting each of the segments in 1 min
                % and calculating the alpha and gammaPower
                alphaPos = intersect(find(freqVals{1}>=freqList{1}(1)),find(freqVals{1}<=freqList{1}(2)));
                sgPos    = intersect(find(freqVals{1}>=freqList{2}(1)),find(freqVals{1}<=freqList{2}(2)));
                fgPos    = intersect(find(freqVals{1}>=freqList{3}(1)),find(freqVals{1}<=freqList{3}(2)));
                fgPos    = fgPos(fgPos~=(find(freqVals{1}==50)));

                meanPSDValsCurrSeg=cell(1,5);

                if size(psdVals{i},2)==180
                    % to fix a special case where only 180 trails for G1
                    psdVals{i}= psdVals{i}(:,1:120);
                end
                if size(psdVals{i},2)==358
                    % to fix a special case where only 258 trails
                    psdVals{i}(:,359:360)=zeros(size(psdVals{i},1),2);
                    badTrialsList{i}=[badTrialsList{i}; 359; 360];
                end
                startIndList = 1:24:size(psdVals{i},2);
                increment = 23;
                for s=1:length(startIndList)
                    meanPSDValsCurrSeg{s}  = mean(psdVals{i}(:,setdiff(startIndList(s):startIndList(s)+increment,badTrialsList{i})),2);
                    meanAlphaCurrSeg(i,s)  = log10(mean(meanPSDValsCurrSeg{s}(alphaPos),1));
                    meanSGCurrSeg(i,s)     = log10(mean(meanPSDValsCurrSeg{s}(sgPos),1));
                    meanFGCurrSeg(i,s)     = log10(mean(meanPSDValsCurrSeg{s}(fgPos),1));
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    timeVecIni=1;
    dimForTrials = 3; % trials are on the third dimention

    if ~isempty(electrodeList)

        timeRange   = [-1 0];  % baseline range
        protocolName = 'G1';
        protocolIndex = 3;
        [psdValsBl,freqVals{i},psdAcrossElc,tfPowerBl,timeValsTFBL,freqValsTF] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,electrodeList,refScheme,capType,trialAvgFlag,diffTf,badTrials,timeRange);
        numTrials = size(psdValsBl,2);

        if g==6
            % get Baseline Psd:
            meanPSDValsBL  = mean(psdValsBl(:,setdiff(1:numTrials,badTrialsList{protocolIndex})),2);
            alphaBl  = log10(mean(meanPSDValsBL(alphaPos),1));
            SGBl     = log10(mean(meanPSDValsBL(sgPos),1));
            FGBl     = log10(mean(meanPSDValsBL(fgPos),1));
        end

        % get BaseLine TF
        blTfMeanTrials = mean(tfPowerBl(:,:,setdiff(1:numTrials,badTrialsList{1})),dimForTrials); % BL power is calculated from G1 protocol
        blTFMeanTime = mean(blTfMeanTrials,1);
        blTFLog      = conv2Log(blTFMeanTime);
        blTF         = repmat(blTFLog,length(timeValsTFBL),1);

        for i=1:numProtocols
            if i==4||i==8 % for M1 and M2 protocol
                increment=15;
            else
                increment=5;
            end

            % Time-frequency plots
            if ~isempty(psdVals{i})
                numTrials = size(psdVals{i},2);

                if rawTfTrialWise % rawTF plots TF TrialWise for all the trials
                    pcolor(hTF(g,i),1:numTrials,freqVals{i},(psdVals{i}));
                    shading(hTF(g,i),'interp');
                    clim(hTF(g,i),cLimsPlot);
                    ylim(hTF(g,i),freqRangeHz);

                elseif diffTf
                    if any(strcmp(protocolNameList(i),{'EO1','EC1','G1'})) % plot as it is after averaging across trials
                        % plotting index does not change
                        % select only Good trials out of all the trials
                        currSegTFAccTrials = mean(tfPower{i}(:,:,setdiff(1:numTrials,badTrialsList{i})),dimForTrials);
                        pcolor(hTF(g,i),timeValsTF,freqValsTF,10*(log10(currSegTFAccTrials)-blTF)');
                        shading(hTF(g,i),'interp');
                        clim(hTF(g,i),cLimsPlot);
                        ylim(hTF(g,i),freqRangeHz);

                    elseif any(strcmp(protocolNameList(i),{'M1','M2'})) % change the index and plot segment wise. TA
                        % creating new indexes for M1 and M1 (a,b and C)
                        if any(strcmp(protocolNameList(i),{'M1'}))
                            j=i;
                            modIndexForMedSeg = j;
                        else % for M2
                            j=10;
                            modIndexForMedSeg = j;
                        end
                        startTrialInd=1;
                        trialIncrement=119;

                        % for a speacial case where the numtrials is less than 360 for meditation segemnts
                        if size(tfPower{i},3)<360
                            tfPower{i}(:,:,359:360)=repmat(zeros(size(tfPower{i},1),size(tfPower{i},2)),[1 1 2]);
                            badTrialsList{i}=[badTrialsList{i}; 359; 360];
                        end

                        for medSegIndex=1:3
                            endTrialIndex = startTrialInd+trialIncrement;
                            currSegTFAllTrials = tfPower{i}(:,:,setdiff(startTrialInd:endTrialIndex,badTrialsList{i}));
                            currSegTFAccTrials = mean(currSegTFAllTrials,dimForTrials);

                            pcolor(hTF(g,modIndexForMedSeg),timeValsTF,freqValsTF,10*(log10(currSegTFAccTrials)-blTF)');
                            shading(hTF(g,modIndexForMedSeg),'interp');
                            clim(hTF(g,modIndexForMedSeg),cLimsPlot);
                            ylim(hTF(g,modIndexForMedSeg),freqRangeHz);

                            modIndexForMedSeg = j + medSegIndex;
                            startTrialInd     = startTrialInd+120;

                            if g==1 && medSegIndex==3  && any(strcmp(protocolNameList(i),{'M2'}))
                                colorBarAx=colorbar(hTF(g,modIndexForMedSeg-1));
                                colorBarAx.Position(1)=colorBarAx.Position(1)+0.025;
                            end
                        end

                    elseif any(strcmp(protocolNameList(i),{'G2','EO2','EC2'}))
                        modIndexForAFterMedSeg = i+2; % new Indexes after M1 (shifted by 2)
                        % select only Good trials out of all the trials
                        currSegTFAccTrials = mean(tfPower{i}(:,:,setdiff(1:size(tfPower{1},3),badTrialsList{i})),dimForTrials);
                        pcolor(hTF(g,modIndexForAFterMedSeg),timeValsTF,freqValsTF,10*(log10((currSegTFAccTrials))-blTF)');
                        shading(hTF(g,modIndexForAFterMedSeg),'interp');
                        clim(hTF(g,modIndexForAFterMedSeg),cLimsPlot);
                        ylim(hTF(g,modIndexForAFterMedSeg),freqRangeHz);
                    end

                else % plots diffTF for all the trials (no badTrial rejection)
                    bl = repmat(meanPSDVals{1},1,numTrials);
                    pcolor(hTF(g,i),1:numTrials,freqVals{i},10*(log10((psdVals{i}))-log10(bl)));
                    shading(hTF(g,i),'interp');
                    clim(hTF(g,i),cLimsPlot);
                    ylim(hTF(g,i),freqRangeHz);
                end

                hold(hTF(g,i),'on');

                if ~isempty(badTrialsList{i})
                    % showing badTrials and total elecs in eachGroup
                    if rawTfTrialWise
                        plot(hTF(g,i),badTrialsList{i},freqRangeHz(2)-1,'k.');
                    end
                    text(1,freqRangeHz(2)-5,['N=' num2str(numGoodElectrodesList(i))],'parent',hTF(g,i));
                end
            end

            if (i==1 && g<numGroups)
                set(hTF(g,i),'XTickLabel',[]); % only remove x label
            elseif (i>1 && g<numGroups)
                set(hTF(g,i),'XTickLabel',[],'YTickLabel',[]);
            elseif (i>1 && g==numGroups)
                set(hTF(g,i),'YTickLabel',[]); % only remove y label
            end

            if g==6
                % plot the power vs time for alpha,sg and fg
                timeVec=timeVecIni:(timeVecIni+(increment-1));
                plot(hPowerVsTime(1),timeVec,10*(nonzeros(meanAlphaCurrSeg(i,:))-repmat(alphaBl,length(nonzeros(meanAlphaCurrSeg(i,:))),1)),'-o','MarkerSize',6,'MarkerFaceColor',colorNames{i},'MarkerEdgeColor',colorNames{i},'color',colorNames{i});
                set(hPowerVsTime(1),'XTickLabel',[]); hold(hPowerVsTime(1),'on');

                plot(hPowerVsTime(2),timeVec,10*(nonzeros(meanSGCurrSeg(i,:))-repmat(SGBl,length(nonzeros(meanSGCurrSeg(i,:))),1)),'-o','MarkerSize',6,'MarkerFaceColor',colorNames{i},'MarkerEdgeColor',colorNames{i},'color',colorNames{i});
                set(hPowerVsTime(2),'XTickLabel',[]); hold(hPowerVsTime(2),'on');

                plot(hPowerVsTime(3),timeVec,10*(nonzeros(meanFGCurrSeg(i,:))-repmat(FGBl,length(nonzeros(meanFGCurrSeg(i,:))),1)),'-o','MarkerSize',6,'MarkerFaceColor',colorNames{i},'MarkerEdgeColor',colorNames{i},'color',colorNames{i});
            end
            timeVecIni=timeVecIni+increment;
        end

        ylabel(hTF(g,1),groupNameList{g});
        xlabel(hPowerVsTime(3),'Time (min)');
        ylabel(hPowerVsTime(3),'\Delta Power (dB)'); hold(hPowerVsTime(3),'on');

        % meanPSD Plots
        for i = 1:numPSDComparisons
            for s=1:length(comparePSDConditions{i})
                conditionNum = comparePSDConditions{i}(s);
                if rawTfTrialWise
                    plot(hPSD(g,i),freqVals{conditionNum},meanPSDVals{conditionNum},'color',colorNames{conditionNum});
                else
                    if ~isempty(meanPSDVals{conditionNum})
                        % bl = repmat(meanPSDVals{1},1,numTrials); % including all the trials
                        plot(hPSD(g,i),freqVals{conditionNum},10*(log10(meanPSDVals{conditionNum})-log10(meanPSDVals{1})),'color',colorNames{conditionNum});
                    end
                end
                hold(hPSD(g,i),'on');
            end

            if ~rawTfTrialWise
                plot(hPSD(g,i),freqVals{conditionNum},zeros(1,length(freqVals{conditionNum})),'k--');
            end

            xlim(hPSD(g,i),freqRangeHz);
            ylim(hPSD(g,i),cLimsPlot);
            if g<numGroups
                set(hPSD(g,i),'XTickLabel',[]);
            end

            if g==1
                title(hPSD(g,i),comparePSDConditionStr{i});
            end
        end
    end % if
end

if rawTfTrialWise
    xlabel(hTF(numGroups,1),'Trials');
else
    xlabel(hTF(numGroups,1),'Time(s)');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting data for the Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psdAcrossElcAccProtocol = cell(1,numProtocols);
timeRange   = [-1 0];  %baseline range
protocolName = 'G1';
[~,freqVals{i},psdAcrossElcBL,~,~,~] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,allElectrodeList,refScheme,capType,trialAvgFlag,diffTf,badTrials,timeRange);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    [~,freqVals,psdAcrossElc] = getData(subjectName,expDate,protocolName,folderSourceString,gridType,allElectrodeList,refScheme,capType,trialAvgFlag,diffTf,badTrials);
    psdAcrossElcAccProtocol{i}  = psdAcrossElc;
end

numElectrodes = length(allElectrodeList);
alphaPos = intersect(find(freqVals>=freqList{1}(1)),find(freqVals<=freqList{1}(2)));
sgPos    = intersect(find(freqVals>=freqList{2}(1)),find(freqVals<=freqList{2}(2)));
fgPos    = intersect(find(freqVals>=freqList{3}(1)),find(freqVals<=freqList{3}(2)));
fgPos    = fgPos(fgPos~=(find(freqVals==50)));

for i=1:numProtocols
    psdDataThisProtocol = psdAcrossElcAccProtocol{i};
    chMeanAlphaThisElec=[];
    chMeanSGThisElec =[];
    chMeanFGThisElec = [];
    noseDir = '+X';

    for e=1:numElectrodes
        psdDataThisElec = squeeze(psdDataThisProtocol(e,:,:));
        psdDataThisElecGoodTrials = mean(psdDataThisElec(:,setdiff(1:size(psdDataThisElec,2),badTrialsList{i})),2);
        psdDataThisElecBL = squeeze(psdAcrossElcBL(e,:,:));
        psdDataThisElecGoodTrialsBL = mean(psdDataThisElecBL(:,setdiff(1:size(psdDataThisElecBL,2),badTrialsList{i})),2);

        % get BLdata:
        BLmeanAlphaThisElec   = log10(mean(psdDataThisElecGoodTrialsBL(alphaPos),1));
        BLmeanSGThisElec      = log10(mean(psdDataThisElecGoodTrialsBL(sgPos),1));
        BLmeanFGThisElec      = log10(mean(psdDataThisElecGoodTrialsBL(fgPos),1));

        % get RawData:
        meanAlphaThisElec   = log10(mean(psdDataThisElecGoodTrials(alphaPos),1));
        meanSGThisElec      = log10(mean(psdDataThisElecGoodTrials(sgPos),1));
        meanFGThisElec      = log10(mean(psdDataThisElecGoodTrials(fgPos),1));

        % change in Power:
        chMeanAlphaThisElec(e,1) = 10*(meanAlphaThisElec-BLmeanAlphaThisElec) ;
        chMeanSGThisElec(e,1)    = 10*(meanSGThisElec- BLmeanSGThisElec) ;
        chMeanFGThisElec(e,1)    = 10*(meanFGThisElec - BLmeanFGThisElec) ;
    end

    % plot topoplot
    x = load([capType 'Labels.mat']);
    montageLabels = x.montageLabels(:,2);
    x = load([capType '.mat']);
    montageChanlocs = x.chanlocs;

    if i>4
        j=i+2;
        set(hTopo(1:3,5),'visible','off');
        set(hTopo(1:3,6),'visible','off');
        set(hTopo(1:3,11),'visible','off');
        set(hTopo(1:3,12),'visible','off');
    elseif i==1
        j=i;
        text(-1.5,0,'\Delta alpha','Parent', hTopo(1,1));
        text(-1.5,0,'\Delta SG','Parent', hTopo(2,1));
        text(-1.5,0,'\Delta FG','Parent', hTopo(3,1));
    else
        j=i;
    end
    axes(hTopo(1,j)); topoplot(chMeanAlphaThisElec,montageChanlocs,'electrodes','off'); clim([-5 5]);
    axes(hTopo(2,j)); topoplot(chMeanSGThisElec,montageChanlocs,'electrodes','off');    clim([-5 5]);
    axes(hTopo(3,j)); topoplot(chMeanFGThisElec,montageChanlocs, 'electrodes','off');    clim([-5 5]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Function End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psd,freqVals,psdAllElc,tfPower,timeValsTF,freqValsTF] = getData(subjectName,expDate,protocolName, ...
    folderSourceString,gridType,electrodeList,refScheme,capType,trialAvgFlag,tfFlag,badTrials,timeRange)
% gives the powerData for all the trials and electrodes

if ~exist('TFFlag','var')    || isempty(tfFlag); tfFlag= 1; end
if ~exist('badTrials','var') || isempty(tfFlag); badTrials= []; end
if ~exist('timeRange','var') || isempty(timeRange); timeRange= [0.25 1.25]; end

% timeRange   = [-0.25 1.25];
tapers      = [1 1];
freqRange   = [0 100];

disp(protocolName);

folderExtract = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'extractedData');
folderSegment = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData');

if ~exist(folderExtract,'file')
    disp([folderExtract ' does not exist']);
    psd = []; freqVals=[];
else
    numElectrodes = length(electrodeList);

    t = load(fullfile(folderSegment,'LFP','lfpInfo.mat'));
    timeVals = t.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    goodTimePos = find(timeVals>=timeRange(1),1) + (1:round(Fs*diff(timeRange)));

    % Setting up multitaper
    params.tapers   = tapers;
    params.pad      = -1; % no padding
    params.Fs       = Fs;
    params.fpass    = freqRange;
    winSize         = 0.25;
    winStep         = 0.025; % 4Hz resolution
    movingwin       = [winSize winStep];

    if trialAvgFlag
        params.trialave = 1;
    else
        params.trialave = 0;
    end

    for i=1:numElectrodes
        if refScheme==2
            load(['bipChInfo' capType '.mat']); % loads the bipolar list
            analogElecs = bipolarLocs(electrodeList(i),:);
            e1 = load(fullfile(folderSegment,'LFP',['elec' num2str(analogElecs(1)) '.mat']));
            e2 = load(fullfile(folderSegment,'LFP',['elec' num2str(analogElecs(2)) '.mat']));
            e.analogData = e1.analogData-e2.analogData; % updated analogData for the bipolar
        else % default: uniPolar
            e = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeList(i)) '.mat']));
        end
        params.trialave = 0;
        [psdTMP(i,:,:),freqVals] = mtspectrumc(e.analogData(:,goodTimePos)',params); %#ok<AGROW>
        [tfPowerTMP(i,:,:,:),timeValsTF0,freqValsTF] = mtspecgramc(e.analogData(:,goodTimePos)',movingwin,params); % for all the trials
        timeValsTF = timeValsTF0 + timeRange(1);
    end
    psd = squeeze(mean(psdTMP,1));          % mean across Electrode psd
    psdAllElc = psdTMP;
    tfPower = squeeze(mean(tfPowerTMP,1));  % mean across electrodes. Already the good electrodes were selected
end
end

function displayBadElecs(hBadElectrodes,subjectName,expDate,protocolName,folderSourceString,gridType,capType,badTrialNameStr,badElectrodes,hStats)

if ~exist('gridType','var');        gridType = 'EEG';                   end
if ~exist('capType','var');         capType = 'actiCap64_2019';         end
if ~exist('badTrialNameStr','var'); badTrialNameStr = '_v5';            end
if ~exist('badElectrodes','var');   badElectrodes = [];                 end
if ~exist('hStats','var');          hStats = [];                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folderSegment = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName,'segmentedData');
badTrialsInfo = load(fullfile(folderSegment,['badTrials' badTrialNameStr '.mat']));

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']); montageChanlocs = x.chanlocs;

if ~isfield(badTrialsInfo,'eegElectrodeLabels') % Check if the labels match with the save labels, if these labels have been saved
    disp('Electrode labels not specified in badTrials file. Taking from Montage...');
else
    if ~isequal(montageLabels(:),badTrialsInfo.eegElectrodeLabels(:))
        error('Montage labels do not match with channel labels in badTrials');
    else
        disp('Montage labels match with channel labels in badTrials');
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
declaredBadElerodeColor = 'k';
if isempty(badElectrodes)
    badImpedanceElectrodes = badTrialsInfo.badElecs.badImpedanceElecs;
    noisyElectrodes = badTrialsInfo.badElecs.noisyElecs;
    flatPSDElectrodes = badTrialsInfo.badElecs.flatPSDElecs;
    declaredBadElecs  =  getDeclaredBadElecs(0.35);
else
    badImpedanceElectrodes = badElectrodes.badImpedanceElecs;
    noisyElectrodes = badElectrodes.noisyElecs;
    flatPSDElectrodes = badElectrodes.flatPSDElecs;
    declaredBadElecs = badElectrodes.declaredBadElecs ;
end

topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{badImpedanceElectrodes,'o',badImpedanceElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{noisyElectrodes,'o',noisyElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{flatPSDElectrodes,'o',flatPSDElectrodeColor,electrodeSize});
topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{declaredBadElecs,'o',declaredBadElerodeColor,electrodeSize});
%topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{highPriorityElectrodeList,'o',highPriorityElectrodeColor,electrodeSize});
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

if ~isempty(hStats)
    axes(hStats)
    set(hStats,'visible','off');
    text(0.05,0.9,'bad Impedance','color',badImpedanceElectrodeColor);
    text(0.05,0.75,num2str(badImpedanceElectrodes(:)'),'color',badImpedanceElectrodeColor);

    text(0.05,0.6,'Noisy','color',noisyElectrodeColor);
    text(0.05,0.45,num2str(noisyElectrodes(:)'),'color',noisyElectrodeColor);

    text(0.05,0.3,'FlatPSD','color',flatPSDElectrodeColor);
    text(0.05,0.15,num2str(flatPSDElectrodes(:)'),'color',flatPSDElectrodeColor);

end
end

function allBadElecs = getAllBadElecs(badElectrodes)
allBadElecs = [badElectrodes.badImpedanceElecs; badElectrodes.noisyElecs; badElectrodes.flatPSDElecs];
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
    topoplot(zeros(1,numElectrodes),montageChanlocs,'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize});
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank');

axes(hPlots(2))
set(hPlots(2),'visible','off');
for i=1:numElectrodeGroups
    text(0.05,0.9-0.15*(i-1),groupNameList{i},'color',electrodeGroupColorList(i,:),'unit','normalized');
end
end