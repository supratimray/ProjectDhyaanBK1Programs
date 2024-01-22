% display Time-frequency data individually for all subjects
clear; clf;
comparisonStr = 'paired';
protocolNameList = [{'G1'} {'G2'} {'M2'}];

badEyeCondition = 'ep';
badTrialVersion = 'v8';
cutoffNumElectrodes = 3;

freqRange = [24 34]; % Show this range
timeLims = [-0.5 1];
freqLims = [0 70];
cLims = [-5 5];

diffTF = 1;
badElecRejectionFlag = 1;
baselineRange = [-1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get subject list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(comparisonStr,'paired')
    pairedSubjectNameList = getPairedSubjectsBK1;
    subjectNameListsTMP{1} = pairedSubjectNameList(:,1);
    subjectNameListsTMP{2} = pairedSubjectNameList(:,2);
else
    [~, meditatorList, controlList] = getGoodSubjectsBK1;
    subjectNameListsTMP{1} = meditatorList;
    subjectNameListsTMP{2} = controlList;
end

%%%%%%%%%%%%%%%%%% Next, divide into males and females %%%%%%%%%%%%%%%%%%%%
[subjectNameListAll,~,~,ageListAll,genderListAll] = getDemographicDetails('BK1');

subjectNameListGender{1} = subjectNameListAll(strcmpi(genderListAll, 'M'));
subjectNameListGender{2} = subjectNameListAll(strcmpi(genderListAll, 'F'));

subjectNameLists = cell(2,2);
numSubjects = zeros(2,2);
ageList = cell(2,2);
for i=1:2 % Gender
    for j=1:2 % Meditator/Control
        % Get subjects sorted by age
        tmpList = intersect(subjectNameListGender{i},subjectNameListsTMP{j});
        numSubjects(i,j) = length(tmpList);

        tmpAges = zeros(1,numSubjects(i,j));
        for k=1:numSubjects(i,j)
            tmpAges(k) = ageListAll(strcmp(tmpList{k},subjectNameListAll));
        end
        [~,pos] = sort(tmpAges);
        subjectNameLists{i,j} = tmpList(pos);
        ageList{i,j} = tmpAges(pos);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numProtocols = length(protocolNameList);
maxNumSubjects = max(numSubjects(:));
hTF = getPlotHandles(maxNumSubjects,4*numProtocols,[0.05 0.05 0.9 0.9],0.005,0.005);

genderStr = [{'M'} {'F'}];
subjectTypeStr = [{'Med'} {'Con'}];

for i=1:2 % Gender
    for j=1:2 % Meditator/Control
        columnIndexList = 2*numProtocols*(i-1) + numProtocols*(j-1) + (1:numProtocols);
        subjectNamesTMP = subjectNameLists{i,j};
        numSubjectsThisCondition = length(subjectNamesTMP);

        for k=1:numSubjectsThisCondition
            subjectName = subjectNamesTMP{k};
            disp(subjectName);
            fileName = fullfile('savedData',[subjectName '_' badEyeCondition '_' badTrialVersion '_TF.mat']);
            tmpData = load(fileName);

            for p=1:numProtocols
                hPlot = hTF(k,columnIndexList(p));

                badElectrodes = tmpData.badElectrodes{p};
                goodElecsInd    = not(ismember(tmpData.electrodeList,badElectrodes));
                
                if sum(goodElecsInd)>=cutoffNumElectrodes
                    meanTFPower = squeeze(mean(tmpData.tfPower{p}(goodElecsInd,:,:),1));
                    logP = log10(meanTFPower);
                    baselinePower = mean(logP(tmpData.timeValsTF>=baselineRange(1) & tmpData.timeValsTF<=baselineRange(2),:));
                    
                    if diffTF
                        pcolor(hPlot,tmpData.timeValsTF,tmpData.freqValsTF,10*(logP'- repmat(baselinePower',1,length(tmpData.timeValsTF))));
                    else
                        pcolor(hPlot,tmpData.timeValsTF,tmpData.freqValsTF,logP'); %#ok<UNRCH>
                    end

                    % Add labels
                    shading(hPlot,'interp');
                    line(timeLims,[freqRange(1) freqRange(1)],'color','k','parent',hPlot);
                    line(timeLims,[freqRange(2) freqRange(2)],'color','k','parent',hPlot);

                    axis(hPlot,[timeLims freqLims]);
                    clim(hPlot,cLims);
                    if p==1
                        text(timeLims(1),freqLims(2)-10,[subjectName ',N=' num2str(tmpData.numTrials(p))],'Parent',hPlot);
                    elseif p==2
                        text(timeLims(1),freqLims(2)-10,[num2str(ageList{i,j}(k)) 'Yr, N=' num2str(tmpData.numTrials(p))],'Parent',hPlot);
                    else
                        text(timeLims(1),freqLims(2)-10,['N=' num2str(tmpData.numTrials(p))],'Parent',hPlot);
                    end
                    set(hPlot, 'TickDir', 'out');
                    set(hPlot,'fontsize',10);
                    set(hPlot,'fontweight','bold');

                    if k==numSubjectsThisCondition
                        xlabel(hPlot,'Time(s)');
                    else
                        set(hPlot,'Xticklabel',[]);
                    end

                    if columnIndexList(p)==1
                        if k==numSubjectsThisCondition
                            ylabel(hPlot,'Frequency (Hz)');
                        end
                    else
                        set(hPlot,'Yticklabel',[]);
                    end

                    if k==1
                        title(hPlot,[genderStr{i} ', ' subjectTypeStr{j} ',' protocolNameList{p}]);
                    end
                else
                    disp(['Not enough good electrodes to plot for ' subjectName tmpData.protocolNameList{p}]);
                end
            end
        end

        % Hide remaining plots
        set(hTF(numSubjectsThisCondition+1:maxNumSubjects,columnIndexList), 'visible', 'off');
    end
end

% set the colorbar
colormap jet
axes(hTF(1,1));
hc = colorbar('Position', [0.96 0.05 0.01 0.1]);
hc.FontSize         = 10;
hc.Label.FontSize   = 10;
hc.Label.FontWeight = 'bold';
hc.Label.String = ['\Delta Power' '(dB)'];