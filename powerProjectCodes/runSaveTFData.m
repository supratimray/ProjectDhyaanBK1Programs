% runSavePowerData
% This program computes the TF data for baseline and stimulus periods for
% all the protocols for all the subjects and saves them locally.

folderSourceString = 'N:\Projects\ProjectDhyaan\BK1'; % Segmented data which is used for all analysis is kept at {folderSourceString}\data\segmentedData

goodSubjectList = getGoodSubjectsBK1; % subjects for which data is computed
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

saveFolderName = 'savedData'; % data is saved locally in this folder

%%%%%%%%%%%%%%% Choices for trial rejection based on eye data %%%%%%%%%%%%%
% The program used to find bad electrodes and a common set of bad trials is
% first run without considering any information about eye position which is
% available through the eye tracker. This is the first level analysis, for
% which the data is saved as badTrials_wo_{badTrialVersion}

% Subsequently, eye data is analyzed and trials in which eye position is
% out of a prespecified window are saved as badTrials_ep

% Finally, eye movements (microsaccades) are analyzed and trials with
% microsaccades are saved as badTrials_em

% If {badEyeCondition} is chosen as 'ep' or 'em', we first find bad trials
% from the 'wo' condition and then appends the extra bad trials based on
% eye data analysis.

badEyeCondition = 'ep'; % use 'wo' when eye data is not used, 'ep' for eye position and 'em' for eye movement
badTrialVersion = 'v8';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
useTheseIndices = 1:length(goodSubjectList);
timeRange = [-1 1.25];

%%%%%%%%%%%%%%%%%%%%%% Get TF data for the protocols where stimulus was presented %%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';
elecRegion = 'highPriority';
protocolNameList = [{'G1'} {'G2'} {'M2'}];
refScheme = 1;   % for unipolar
[~,~,~,electrodeGroupList,groupNameList,highPriorityElectrodeNums] = electrodePositionOnGrid(1,gridType,[],capType,refScheme);
numGroups = length(electrodeGroupList);
electrodeGroupList{numGroups+1} = highPriorityElectrodeNums;
groupNameList{numGroups+1} = 'highPriority';
numGroups = numGroups+1;
electrodeList = electrodeGroupList{strcmp(groupNameList,elecRegion)};

for i=1:length(useTheseIndices)
    subjectName = goodSubjectList{useTheseIndices(i)};
    disp(['Analyzing for subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};

    [tfPower,timeValsTF,freqValsTF,numTrials,badElectrodes] = getTFData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,timeRange,electrodeList);

    makeDirectory(saveFolderName);
    fileNameSave = fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' 'TF']);
    save(fileNameSave,'tfPower','timeValsTF','freqValsTF','numTrials','electrodeList','badElectrodes','timeRange','protocolNameList');
end