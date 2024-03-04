% This program saves data in fieldtrip format

function saveFTData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,ftDataFolder)

numProtocols = length(protocolNameList);

for i=1:numProtocols
    protocolName = protocolNameList{i};
    saveFTDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion,ftDataFolder);
end
end

function saveFTDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion,ftDataFolder)

gridType = 'EEG'; capType = 'actiCap64_UOL';

% Fieldtrip folder
ftDataFolder = fullfile(ftDataFolder,subjectName);
makeDirectory(ftDataFolder);
ftDataFileName = fullfile(ftDataFolder,[protocolName '_' badEyeCondition '_' badTrialVersion]);

timeValsFileName = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData','LFP','lfpInfo.mat');
if ~exist(timeValsFileName,'file')
    data = []; numGoodTrials = 0;
else
    x = load(timeValsFileName);
    timeVals = x.timeVals;

    % Cap details
    elec = getElectrodeDetails_ft(capType);
    numElectrodes = length(elec.label);

    % Get bad electrodes and trials
    [badTrials,badElectrodes] = getBadTrialsAndElectrodes(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion);

    % Get data from all electrodes, including electrodes labelled bad.
    for i=1:numElectrodes
        dataFileName = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData','LFP',['elec' num2str(i) '.mat']);
        x = load(dataFileName);
        sig_all(i,:,:) = x.analogData; %#ok<AGROW>
    end
    sig_all(:,badTrials,:) = []; % Remove bad trials
    numTrials = size(sig_all,2);

    eeg_all=cell(1,numTrials);
    for tr=1:numTrials
        eeg_all{tr}=squeeze(sig_all(:,tr,:));
    end

    chantype=cell(numElectrodes,1);
    chanunit=cell(numElectrodes,1);
    for i=1:numElectrodes
        chantype{i}='eeg';
        chanunit{i}='uV';
    end

    hdr.Fs=round(1/(timeVals(2)-timeVals(1)));
    hdr.nChans=numElectrodes;
    hdr.nSamples=length(timeVals);
    hdr.nTrials=numTrials;
    hdr.label=elec.label;
    hdr.chantype=chantype;
    hdr.chanunit=chanunit;

    times=cell(1,numTrials);
    for i=1:length(times)
        times{i}=timeVals;
    end

    data.hdr=hdr;
    data.label=hdr.label;
    data.trial=eeg_all; % cell of each trial as (chans x timpts)
    data.time=times; % cell of each trial timepoints
    data.fsample=hdr.Fs;
    data.cfg=[];
    data.elec=elec;

    data.timeVals = timeVals;
    data.numTrials = numTrials;
    data.badElecs = badElectrodes;

    numGoodTrials = data.numTrials;
end
save(ftDataFileName,'data','numGoodTrials');
end
function elec = getElectrodeDetails_ft(capType)

if strcmp(capType,'actiCap64_UOL')

    labels = load([capType 'Labels.mat']);
    chanlocs = load([capType '.mat']); %#ok<NASGU>

    dimLabels = 'XYZ';
    for dim = 1:3
        for e = 1:64
            clocs(e,dim) = eval(['chanlocs.chanlocs(' int2str(e) ').' dimLabels(dim)]); %#ok<AGROW>
        end
    end
    elec.chanpos = clocs;
    elec.chantype = repmat({'eeg'},64,1);
    elec.elecpos = clocs;
    elec.label = labels.montageLabels(:,2);
    elec.type = 'eeg1010';
    elec.unit = 'cm';
else
    disp([capType{1} ' not specified']);
    elec = [];
end
end