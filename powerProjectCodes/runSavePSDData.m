% runSavePowerData
% This program computes the PSD data for baseline and stimulus periods for
% all the protocols and all the subjects and saves them locally.

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
stRange = [0.25 1.25];

%%%%%%%%%%%%%%%%%%%%%% Get PSD data for all protocols %%%%%%%%%%%%%%%%%%%%%
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 

for i=1:length(useTheseIndices)
    subjectName = goodSubjectList{useTheseIndices(i)};
    disp(['Analyzing for subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};
    
    [psdValsST,psdValsBL,freqVals,numTrials,badElectrodes] = getPSDData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,stRange);

    makeDirectory(saveFolderName);
    fileNameSave = fullfile(saveFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]);
    save(fileNameSave,'psdValsST','psdValsBL','freqVals','numTrials','badElectrodes');
end