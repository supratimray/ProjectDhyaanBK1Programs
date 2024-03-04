function [badTrials,badElectrodes] = getBadTrialsAndElectrodes(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion)
gridType = 'EEG';

badTrialNameStr = ['_wo_' badTrialVersion];
badFileName = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData',['badTrials' badTrialNameStr '.mat']);
if exist(badFileName,'file')
    x=load(badFileName);
    badTrials = x.badTrials(:);
    badElectrodes = getAllBadElecs(x.badElecs);

    if ~strcmp(badEyeCondition,'wo')
        badFileName2 = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData',['badTrials_' badEyeCondition '.mat']);
        x=load(badFileName2);
        badTrials = unique(cat(1,badTrials,x.badEyeTrials));
    end
else
    disp([badFileName ' does not exist']);
    badTrials=[];
    badElectrodes=[];
end
end

function allBadElecs = getAllBadElecs(badElectrodes)
if exist('getDeclaredBadElecs','file')
    declaredBadElecs = getDeclaredBadElecs';
else
    declaredBadElecs = [];
end

if ~isempty(badElectrodes)
    allBadElecs = unique([badElectrodes.badImpedanceElecs; badElectrodes.noisyElecs; badElectrodes.flatPSDElecs; declaredBadElecs]);
else
    allBadElecs = [];
end
end