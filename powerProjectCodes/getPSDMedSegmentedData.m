function [psdValsST,psdValsBL,freqVals,numTrials,badElectrodes] = getPSDMedSegmentedData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,stRange)
% getPSDMedSegmentedData Compute PSD data for meditation protocols divided into segments
% Processes M1 and M2 protocols in three segments

% Input validation
if ~exist('folderSourceString','var'); folderSourceString=[]; end
if ~exist('badEyeCondition','var'); badEyeCondition = 'wo'; end
if ~exist('badTrialVersion','var'); badTrialVersion = 'v8'; end
if ~exist('stRange','var'); stRange = [0.25 1.25]; end

if isempty(folderSourceString)
    folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numProtocols = length(protocolNameList);
psdValsST = cell(1,numProtocols);
psdValsBL = cell(1,numProtocols);
numTrials = zeros(1,numProtocols);
badElectrodes = cell(1,numProtocols);

% Parameters
trialsPerSegment = 120;

for i=1:numProtocols
    protocolName = protocolNameList{i};
    fprintf('Analyzing protocol segment %s (%d/%d)...\n', protocolName, i, numProtocols);

    % Determine parent protocol and segment
    if contains(protocolName, 'M1')
        parentProtocol = 'M1';
    else
        parentProtocol = 'M2';
    end

    % Get segment index (a=1, b=2, c=3)
    segmentLetter = protocolName(end);
    segmentIndex = double(segmentLetter) - double('a') + 1;

    % Calculate fixed trial ranges for each segment
    startTrial = (segmentIndex-1) * trialsPerSegment + 1;  % 1, 121, or 241
    endTrial = segmentIndex * trialsPerSegment;            % 120, 240, or 360

    fprintf('Processing trials %d to %d for %s\n', startTrial, endTrial, protocolName);

    % Get bad trials from parent protocol
    [badTrials,badElectrodes{i}] = getBadTrialsAndElectrodes(subjectName,expDate,parentProtocol,folderSourceString,badEyeCondition,badTrialVersion);

    % Get PSD data for the segment
    [psdValsST{i},psdValsBL{i},freqVals,numTrials(i)] = getDataForMedSegment(subjectName,expDate,parentProtocol,folderSourceString,stRange,badTrials,startTrial,endTrial);
end
end

function [psdValsST,psdValsBL,freqVals,numTrials] = getDataForMedSegment(subjectName,expDate,protocolName,folderSourceString,stRange,badTrials,startTrial,endTrial)
gridType = 'EEG';
electrodeList = 1:64;
tapers = [1 1];
freqRange = [0 100];
numElectrodes = length(electrodeList);

folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');
timingFile = fullfile(folderSegment,'LFP','lfpInfo.mat');

if ~exist(timingFile,'file')
    disp([timingFile ' does not exist.']);
    psdValsST = []; psdValsBL = []; freqVals = []; numTrials = 0;
    return;
end

% Load timing information
t = load(timingFile);
timeVals = t.timeVals;
Fs = round(1/(timeVals(2)-timeVals(1)));
T = diff(stRange);

% Calculate time indices
goodTimePosST = find(timeVals>=stRange(1),1) + (1:round(Fs*diff(stRange)));
blRange = [-T 0];
goodTimePosBL = find(timeVals>=blRange(1),1) + (1:round(Fs*diff(blRange)));

% Set up multitaper parameters
params.tapers = tapers;
params.pad = -1;
params.Fs = Fs;
params.fpass = freqRange;
params.trialave = 1;

% Load first electrode to get data dimensions
e = load(fullfile(folderSegment,'LFP','elec1.mat'));

% Get segment trials that are good
segmentTrials = startTrial:endTrial;  % Use exact trial range
goodTrials = intersect(segmentTrials, setdiff(1:size(e.analogData,1), badTrials));
numTrials = length(goodTrials);

fprintf('Found %d good trials in range %d-%d\n', numTrials, startTrial, endTrial);

if numTrials > 0
    psdValsST = zeros(numElectrodes, length(params.fpass(1):1/T:params.fpass(2)));
    psdValsBL = zeros(numElectrodes, length(params.fpass(1):1/T:params.fpass(2)));

    for i=1:numElectrodes
        e = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeList(i)) '.mat']));
        [psdValsST(i,:),freqVals] = mtspectrumc(e.analogData(goodTrials,goodTimePosST)',params);
        psdValsBL(i,:) = mtspectrumc(e.analogData(goodTrials,goodTimePosBL)',params);
    end
else
    psdValsST = [];
    psdValsBL = [];
    freqVals = [];
end
end