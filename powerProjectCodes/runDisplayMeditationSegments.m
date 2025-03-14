function runDisplayMeditationSegments

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.125;
colormap jet

%%%%%%%%%%%%%%%%%%%%%%%%%%% Subject Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel1 = uipanel('Title','Subjects','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.025 1-panelHeight 0.15 panelHeight]);

% Comparison - paired or unpaired
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','Comparison','FontSize',fontSizeSmall);
comparisonList = [{'paired'} {'unpaired'}];
hComparison = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',comparisonList,'FontSize',fontSizeSmall);

% Gender - all, male, female
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','Gender','FontSize',fontSizeSmall);
genderList = [{'all'} {'male'} {'female'}];
hGender = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',genderList,'FontSize',fontSizeSmall);

% Age - all, young, mid
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','Age','FontSize',fontSizeSmall);
ageList = [{'all'} {'young'} {'mid'}];
hAge = uicontrol('Parent',hPanel1,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',ageList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%% Protocol Details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel2 = uipanel('Title','Protocol','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.175 1-panelHeight 0.15 panelHeight]);

% Protocol
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0.75 0.5 0.25],'Style','text','String','Protocol','FontSize',fontSizeSmall);
protocolList = [{'M1'} {'M2'}];
hProtocol = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0.75 0.5 0.25],'Style','popup','String',protocolList,'FontSize',fontSizeSmall);

% Analysis Choice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0.5 0.5 0.25],'Style','text','String','Analysis','FontSize',fontSizeSmall);
analysisChoiceList1 = [{'spontaneous (bl)'} {'stimulus (st)'} {'combined'}];
analysisChoiceList2 = [{'bl'} {'st'} {'combined'}];
hAnalysisType = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0.5 0.5 0.25],'Style','popup','String',analysisChoiceList1,'FontSize',fontSizeSmall);

% RefChoice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0.25 0.5 0.25],'Style','text','String','Ref Choice','FontSize',fontSizeSmall);
refChoiceList = [{'none'} {'EO1'} {'EC1'} {'G1'} {'G2'} {'EO2'} {'EC2'} ...
    {'M1a'} {'M1b'} {'M1c'} {'M1'} ...
    {'M2a'} {'M2b'} {'M2c'} {'M2'}];
hRefChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0.25 0.5 0.25],'Style','popup','String',refChoiceList,'FontSize',fontSizeSmall);

% Elec Choice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0 0.5 0.25],'Style','text','String','Elec Choice','FontSize',fontSizeSmall);
elecChoiceList = [{'Occipital'} {'Frontal'}];
hElecChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 0.25],'Style','popup','String',elecChoiceList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Bad Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel3 = uipanel('Title','Bad Elec Condition','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.325 1-panelHeight 0.15 panelHeight]);

% Bad Eye condition
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','BadEyeCondition','FontSize',fontSizeSmall);
badEyeConditionList1 = [{'eye position (ep)'} {'none (wo)'}]; badEyeConditionList2 = [{'ep'} {'wo'}];
hBadEye = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',badEyeConditionList1,'Value',1,'FontSize',fontSizeSmall);

% Bad Trial Version
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','BadTrialVersion','FontSize',fontSizeSmall);
badTrialVersionList = {'v8'};
hBadTrialVersion = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',badTrialVersionList,'FontSize',fontSizeSmall);

% Bad Electrode Choice
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','BadElecChoice','FontSize',fontSizeSmall);
badElectrodeChoiceList = [{'Reject badElectrodes of protocolName'} {'Reject common badElectrodes of all protocols'} {'Reject badElectrodes of G1'}];
hBadElectrodeChoice = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',badElectrodeChoiceList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%% Freq Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel4 = uipanel('Title','Freq Ranges','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.475 1-panelHeight 0.15 panelHeight]);
freqRangeList0{1} = [8 13];
freqRangeList0{2} = [24 34];
freqRangeList0{3} = [35 65];

numFreqRanges = length(freqRangeList0);
hFreqRangeMin = cell(1,numFreqRanges);
hFreqRangeMax = cell(1,numFreqRanges);

for i=1:numFreqRanges
    uicontrol('Parent',hPanel4,'Unit','Normalized','Position',[0 1-i/numFreqRanges 0.5 1/numFreqRanges],'Style','text','String',['Freq Range' num2str(i)],'FontSize',fontSizeSmall);
    hFreqRangeMin{i} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numFreqRanges 0.25 1/numFreqRanges], ...
        'Style','edit','String',num2str(freqRangeList0{i}(1)),'FontSize',fontSizeSmall);
    hFreqRangeMax{i} = uicontrol('Parent',hPanel4,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-i/numFreqRanges 0.25 1/numFreqRanges], ...
        'Style','edit','String',num2str(freqRangeList0{i}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Axis Ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel5 = uipanel('Title','Axis Ranges','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.625 1-panelHeight 0.15 panelHeight]);
axisRangeList0{1} = [0 80]; axisRangeName{1} = 'Freq Lims (Hz)';
axisRangeList0{2} = [-2.5 2.5]; axisRangeName{2} = 'YLims';
axisRangeList0{3} = [-1.5 1.5]; axisRangeName{3} = 'cLims (topo)';

numAxisRanges = length(axisRangeList0);
hAxisRangeMin = cell(1,numAxisRanges);
hAxisRangeMax = cell(1,numAxisRanges);

for i=1:numAxisRanges
    uicontrol('Parent',hPanel5,'Unit','Normalized','Position',[0 1-i/numAxisRanges 0.5 1/numAxisRanges],'Style','text','String',axisRangeName{i},'FontSize',fontSizeSmall);
    hAxisRangeMin{i} = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numAxisRanges 0.25 1/numAxisRanges], ...
        'Style','edit','String',num2str(axisRangeList0{i}(1)),'FontSize',fontSizeSmall);
    hAxisRangeMax{i} = uicontrol('Parent',hPanel5,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.75 1-i/numAxisRanges 0.25 1/numAxisRanges], ...
        'Style','edit','String',num2str(axisRangeList0{i}(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Cutoff Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel6 = uipanel('Title','Cutoffs','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.775 1-panelHeight 0.1 panelHeight]);
cutoffList0 = [3 10 30]; cutoffNames = [{'Num Elecs'} {'Num Trials'} {'Time (min)'}];

numCutoffRanges = length(cutoffList0);
hCutoffs = cell(1,numCutoffRanges);

for i=1:numCutoffRanges
    uicontrol('Parent',hPanel6,'Unit','Normalized','Position',[0 1-i/numCutoffRanges 0.5 1/numCutoffRanges],'Style','text','String',cutoffNames{i},'FontSize',fontSizeSmall);
    hCutoffs{i} = uicontrol('Parent',hPanel6,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numCutoffRanges 0.5 1/numCutoffRanges], ...
        'Style','edit','String',num2str(cutoffList0(i)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Plot Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel7 = uipanel('Title','Plot','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.875 1-panelHeight 0.1 panelHeight]);

hUseMedianFlag = uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 3/5 1 1/5],'Style','togglebutton','String','Use Median','FontSize',fontSizeMedium);
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 2/5 0.5 1/5],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0.5 2/5 0.5 1/5],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 0 1 2/5],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define segment list for meditation protocols
segmentList = [{'a'} {'b'} {'c'}];

% Get plot handles using getPlotHandles
hAllPlots.hPSD = getPlotHandles(1, 5, [0.05 0.62 0.9 0.2], 0.02, 0.02, 1);
hAllPlots.hPower = getPlotHandles(3, 5, [0.05 0.06 0.55 0.5], 0.02, 0.02, 0);
hAllPlots.hTopo = getPlotHandles(3, 6, [0.62 0.06 0.36 0.5], 0.001, 0.02, 1);

    function plot_Callback(~,~)

        %%%%%%%%%%%%%%%%%%%%% Get SubjectLists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        comparisonStr=comparisonList{get(hComparison,'val')};
        if strcmp(comparisonStr,'paired')
            pairedSubjectNameList = getPairedSubjectsBK1;
            subjectNameLists{1} = pairedSubjectNameList(:,1);
            subjectNameLists{2} = pairedSubjectNameList(:,2);
            pairedDataFlag      = 1;
        else
            [~, meditatorList, controlList] = getGoodSubjectsBK1;
            subjectNameLists{1} = meditatorList;
            subjectNameLists{2} = controlList;
            pairedDataFlag      = 0;
        end

        % Sub-select Subjects based on Demographics
        [subjectNameList,~,~,ageListAllSub,genderListAllSub] = getDemographicDetails('BK1');
        % sub-select based on Gender
        genderStr=genderList{get(hGender,'Value')};
        maleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'M'));
        femaleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'F'));
        if  strcmp(genderStr,'male')
            subjectNameLists{1} = intersect(subjectNameLists{1},maleSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},maleSubjectNameList,'stable');
        elseif strcmp(genderStr,'female')
            subjectNameLists{1} = intersect(subjectNameLists{1},femaleSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},femaleSubjectNameList,'stable');
        end

        % sub-select based on Age
        ageStr=ageList{get(hAge,'Value')};
        youngSubjectNameList = subjectNameList(ageListAllSub<40);
        midSubjectNameList = subjectNameList(ageListAllSub>=40);

        if  strcmp(ageStr,'young')
            subjectNameLists{1} = intersect(subjectNameLists{1},youngSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},youngSubjectNameList,'stable');
        elseif strcmp(ageStr,'mid')
            subjectNameLists{1} = intersect(subjectNameLists{1},midSubjectNameList,'stable');
            subjectNameLists{2} = intersect(subjectNameLists{2},midSubjectNameList,'stable');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Get GUI parameters
        comparisonStr = comparisonList{get(hComparison,'val')};
        protocolStr = protocolList{get(hProtocol,'val')};
        analysisChoice = analysisChoiceList2{get(hAnalysisType,'val')};
        refChoice = refChoiceList{get(hRefChoice,'val')};
        badEyeCondition = badEyeConditionList2{get(hBadEye,'val')};
        badTrialVersion = badTrialVersionList{get(hBadTrialVersion,'val')};
        badElectrodeRejectionFlag = get(hBadElectrodeChoice,'val');
        useMedianFlag = get(hUseMedianFlag,'val');
        elecChoice = elecChoiceList{get(hElecChoice,'val')};

        % Get frequency ranges from GUI
        freqRangeList = cell(1,numFreqRanges);
        for j=1:numFreqRanges
            freqRangeList{j} = [str2double(get(hFreqRangeMin{j},'String')) str2double(get(hFreqRangeMax{j},'String'))];
        end

        % Get axis ranges from GUI
        axisRangeList = cell(1,numAxisRanges);
        for j=1:numAxisRanges
            axisRangeList{j} = [str2double(get(hAxisRangeMin{j},'String')) str2double(get(hAxisRangeMax{j},'String'))];
        end

        % Get cutoff values from GUI
        cutoffList = zeros(1,numCutoffRanges);
        for j=1:numCutoffRanges
            cutoffList(j) = str2double(get(hCutoffs{j},'String'));
        end

        % main function
        displayPowerDataAllSubjectsMedSegmented(hAllPlots, subjectNameLists, segmentList, comparisonStr, protocolStr, ...
            analysisChoice, refChoice, badEyeCondition, badTrialVersion, badElectrodeRejectionFlag, ...
            freqRangeList, useMedianFlag, axisRangeList, cutoffList, pairedDataFlag, elecChoice);
    end

    function rescale_Callback(~,~)
        axisLims = [str2double(get(hAxisRangeMin{1},'String')) str2double(get(hAxisRangeMax{1},'String')) ...
            str2double(get(hAxisRangeMin{2},'String')) str2double(get(hAxisRangeMax{2},'String'))];
        cLims = [str2double(get(hAxisRangeMin{3},'String')) str2double(get(hAxisRangeMax{3},'String'))];

        % Rescale PSD plot
        axes(hAllPlots.hPSD);
        axis(axisLims);

        % Rescale topoplots
        for p=1:size(hAllPlots.hTopo,1)
            for j=1:size(hAllPlots.hTopo,2)
                clim(hAllPlots.hTopo(p,j),cLims);
            end
        end
    end

    function cla_Callback(~,~)
        % Clear all standard plots
        claGivenPlotHandle(hAllPlots.hPSD);
        claGivenPlotHandle(hAllPlots.hPower);
        claGivenPlotHandle(hAllPlots.hTopo);

        function claGivenPlotHandle(plotHandles)
            % Get size of plot handle array
            [numRows,numCols] = size(plotHandles);

            % Clear each subplot
            for r=1:numRows
                for j=1:numCols
                    if ishandle(plotHandles(r,j))
                        cla(plotHandles(r,j));
                    end
                end
            end
        end
    end
end