% A GUI to choose different options for displaying data

function runDisplayPowerDataAllSubjects

fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;
backgroundColor = 'w'; panelHeight = 0.125;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Subject Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hPanel1 = uipanel('Title','Subjects','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.025 1-panelHeight 0.15 panelHeight]);

% Comparison - paired or unpaired
uicontrol('Parent',hPanel1,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','Comparison','FontSize',fontSizeSmall);
comparisonList = [{'unpaired'} {'paired'}];
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
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','ProtocolName','FontSize',fontSizeSmall);
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];
hProtocol = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',protocolNameList,'FontSize',fontSizeSmall);

% AnalysisChoice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','text','String','Analysis','FontSize',fontSizeSmall);
analysisChoiceList1 = [{'spontaneous (bl)'} {'stimulus (st)'} {'combined'}];
analysisChoiceList2 = [{'bl'} {'st'} {'combined'}];
hAnalysisChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 1/3 0.5 1/3],'Style','popup','String',analysisChoiceList1,'FontSize',fontSizeSmall);

% RefChoice
uicontrol('Parent',hPanel2,'Unit','Normalized','Position',[0 0 0.5 1/3],'Style','text','String','Ref Choice','FontSize',fontSizeSmall);
refChoiceList = [{'none'} protocolNameList];
hRefChoice = uicontrol('Parent',hPanel2,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 0 0.5 1/3],'Style','popup','String',refChoiceList,'FontSize',fontSizeSmall);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Bad Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel3 = uipanel('Title','Bad Electrode Condition','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.325 1-panelHeight 0.15 panelHeight]);

% Bad Eye condition
uicontrol('Parent',hPanel3,'Unit','Normalized','Position',[0 2/3 0.5 1/3],'Style','text','String','BadEyeCondition','FontSize',fontSizeSmall);
badEyeConditionList1 = [{'eye position (ep)'} {'none (wo)'}]; badEyeConditionList2 = [{'ep'} {'wo'}];
hBadEye = uicontrol('Parent',hPanel3,'Unit','Normalized','BackgroundColor', backgroundColor, 'Position', [0.5 2/3 0.5 1/3],'Style','popup','String',badEyeConditionList1,'FontSize',fontSizeSmall);

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
freqRangeList0{2} = [22 34];
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
axisRangeList0{1} = [0 100]; axisRangeName{1} = 'Freq Lims (Hz)';
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
cutoffList0 = [3 50]; cutoffNames = [{'Num Elecs'} {'Num Trials'}];

numCutoffRanges = length(cutoffList0);
hCutoffs = cell(1,numCutoffRanges);

for i=1:numCutoffRanges
    uicontrol('Parent',hPanel6,'Unit','Normalized','Position',[0 1-i/numCutoffRanges 0.5 1/numCutoffRanges],'Style','text','String',cutoffNames{i},'FontSize',fontSizeSmall);
    hCutoffs{i} = uicontrol('Parent',hPanel6,'Unit','Normalized','BackgroundColor', backgroundColor,'Position',[0.5 1-i/numCutoffRanges 0.5 1/numCutoffRanges], ...
    'Style','edit','String',num2str(cutoffList0(i)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Plot Choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hPanel7 = uipanel('Title','Plot','fontSize',fontSizeLarge,'Unit','Normalized','Position',[0.875 1-panelHeight 0.1 panelHeight]);

hUseMedianFlag = uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 2/3 1 1/3],'Style','togglebutton','String','Use Median','FontSize',fontSizeMedium);
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 1/3 0.5 1/3],'Style','pushbutton','String','Rescale','FontSize',fontSizeMedium,'Callback',{@rescale_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0.5 1/3 0.5 1/3],'Style','pushbutton','String','Clear','FontSize',fontSizeMedium,'Callback',{@cla_Callback});
uicontrol('Parent',hPanel7,'Unit','Normalized','Position',[0 0 1 1/3],'Style','pushbutton','String','plot','FontSize',fontSizeMedium,'Callback',{@plot_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodeGroupList = getElectrodeGroups('EEG','actiCap64_UOL');
numGroups = length(electrodeGroupList);
hAllPlots.hPSD  =  getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.3],0.02,0.02,1);
hAllPlots.hPower = getPlotHandles(numFreqRanges,numGroups,[0.05 0.05 0.6 0.45],0.02,0.02,0);
hAllPlots.hTopo0 = getPlotHandles(1,2,[0.675 0.7 0.3 0.15],0.02,0.02,1);
hAllPlots.hTopo1 = getPlotHandles(1,3,[0.675 0.55 0.3 0.13],0.02,0.02,1);
hAllPlots.hTopo2 = getPlotHandles(numFreqRanges,3,[0.675 0.05 0.3 0.45],0.02,0.02,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot_Callback(~,~)

    %%%%%%%%%%%%%%%%%%%%% Get SubjectLists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        comparisonStr=comparisonList{get(hComparison,'val')};

        if strcmp(comparisonStr,'paired')
            pairedSubjectNameList = getPairedSubjectsBK1;
            subjectNameLists{1} = pairedSubjectNameList(:,1);
            subjectNameLists{2} = pairedSubjectNameList(:,2);
        else
            [~, meditatorList, controlList] = getGoodSubjectsBK1;
            subjectNameLists{1} = meditatorList;
            subjectNameLists{2} = controlList;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        protocolName = protocolNameList{get(hProtocol,'val')};
        analysisChoice = analysisChoiceList2{get(hAnalysisChoice,'val')};
        refChoice = refChoiceList{get(hRefChoice,'val')};

        badEyeCondition = badEyeConditionList2{get(hBadEye,'val')};
        badTrialVersion = badTrialVersionList{get(hBadTrialVersion,'val')};
        badElectrodeRejectionFlag = get(hBadElectrodeChoice,'val');
        
        stRange = [0.25 1.25]; % hard coded for now

        freqRangeList = cell(1,numFreqRanges);
        for ii=1:numFreqRanges
            freqRangeList{ii} = [str2double(get(hFreqRangeMin{ii},'String')) str2double(get(hFreqRangeMax{ii},'String'))];
        end

        axisRangeList = cell(1,numAxisRanges);
        for ii=1:numAxisRanges
            axisRangeList{ii} = [str2double(get(hAxisRangeMin{ii},'String')) str2double(get(hAxisRangeMax{ii},'String'))];
        end

        cutoffList = zeros(1,numCutoffRanges);
        for ii=1:numCutoffRanges
            cutoffList(ii) = str2double(get(hCutoffs{ii},'String'));
        end

        useMedianFlag = get(hUseMedianFlag,'val');

        displayPowerDataAllSubjects(subjectNameLists,protocolName,analysisChoice,refChoice,badEyeCondition,badTrialVersion,badElectrodeRejectionFlag,stRange,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots);
    end
    function cla_Callback(~,~)
        claGivenPlotHandle(hAllPlots.hPSD);
        claGivenPlotHandle(hAllPlots.hPower);
        claGivenPlotHandle(hAllPlots.hTopo0);
        claGivenPlotHandle(hAllPlots.hTopo1);
        claGivenPlotHandle(hAllPlots.hTopo2);
        
        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    cla(plotHandles(ii,j));
                end
            end
        end
    end
    function rescale_Callback(~,~)
        axisLims = [str2double(get(hAxisRangeMin{1},'String')) str2double(get(hAxisRangeMax{1},'String')) str2double(get(hAxisRangeMin{2},'String')) str2double(get(hAxisRangeMax{2},'String'))];
        cLims = [str2double(get(hAxisRangeMin{3},'String')) str2double(get(hAxisRangeMax{3},'String'))];

        rescaleGivenPlotHandle(hAllPlots.hPSD,axisLims);
        rescaleZGivenPlotHandle(hAllPlots.hTopo2,cLims);

        function rescaleGivenPlotHandle(plotHandles,axisLims)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    axis(plotHandles(ii,j),axisLims);
                end
            end
        end
        function rescaleZGivenPlotHandle(plotHandles,cLims)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    clim(plotHandles(ii,j),cLims);
                end
            end
        end
    end

end