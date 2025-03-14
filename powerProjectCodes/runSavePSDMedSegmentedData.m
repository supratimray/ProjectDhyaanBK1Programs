% runSavePSDMedSegmentedData
% This program computes the PSD data for baseline and stimulus periods for
% meditation protocols divided into three segments and saves them locally.

% If {badEyeCondition} is chosen as 'ep' or 'em', we first find bad trials
% from the 'wo' condition and then appends the extra bad trials based on
% eye data analysis.

folderSourceString = 'D:\Projects\ProjectDhyaan\BK1';

goodSubjectList = getGoodSubjectsBK1;
[allSubjectNames,expDateList] = getDemographicDetails('BK1');

saveFolderName = 'savedDataMedSegmented';

% Analysis Parameters
badEyeCondition = 'ep';
badTrialVersion = 'v8';
stRange = [0.25 1.25];

% Create save directory if it doesn't exist
makeDirectory(saveFolderName);

% Define meditation protocol segments
protocolNameList = [{'M1a'} {'M1b'} {'M1c'} {'M2a'} {'M2b'} {'M2c'}];
protocolNameMain = {'M1','M2'};

useTheseIndices = 1:length(goodSubjectList);

for i=1:length(useTheseIndices)
    subjectName = goodSubjectList{useTheseIndices(i)};
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};

    fileNameSave = fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2)) '_medseg.mat']);
    fileNameSegmented = fullfile(folderSourceString,'data','segmentedData',subjectName,'EEG',expDate);

    % Check if data file already exists
    if exist(fileNameSave, 'file')
        fprintf('Data file already exists for subject %s, skipping...\n', subjectName);
        continue;
    else
        % Process data if file doesn't exist
        disp(['Analyzing for subject ' subjectName]);
        [psdValsST,psdValsBL,freqVals,numTrials,badElectrodes] = getPSDMedSegmentedData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,stRange);
        protcolLengthTime = saveTimeData(fileNameSegmented,fileNameSave,protocolNameMain,subjectName);
        % Save data
        save(fileNameSave,'psdValsST','psdValsBL','freqVals','numTrials','badElectrodes','protcolLengthTime');
        fprintf('Data saved successfully for subject %s\n', subjectName);
    end
end

fprintf('\nAnalyzed %d subjects\n', length(useTheseIndices));