function [tfPower,timeValsTF,freqVals,numTrials,badElectrodes] = getTFData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,timeRange,electrodeList)

if ~exist('protocolNameList','var')
    protocolNameList = [ {'G1'} {'G2'} {'M2'}];
end
if ~exist('folderSourceString','var');   folderSourceString=[];         end
if ~exist('badEyeCondition','var');      badEyeCondition = 'wo';        end
if ~exist('badTrialVersion','var');      badTrialVersion = 'v8';        end
if ~exist('timeRange','var');            timeRange = [-1 1];         end
if isempty(folderSourceString)
    folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numProtocols = length(protocolNameList);

tfPower = cell(1,numProtocols);
numTrials = zeros(1,numProtocols);

badElectrodes = cell(1,numProtocols);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    [badTrials,badElectrodes{i}] = getBadTrialsAndElectrodes(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion);
    [tfPower{i},timeValsTF,freqVals,numTrials(i)] = getDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,timeRange,badTrials,electrodeList);
end
end

function [tfPower,timeValsTF,freqVals,numTrials] = getDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,timeRange,badTrials,electrodeList)

if ~exist('electrodeList','var');   electrodeList=1:64;         end

gridType = 'EEG';
tapers = [1 1];
freqRange = [0 100];

numElectrodes = length(electrodeList);

folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');
timingFile = fullfile(folderSegment,'LFP','lfpInfo.mat');

if exist(timingFile,'file')

    t = load(timingFile);
    timeVals = t.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    goodTimePos = find(timeVals>=timeRange(1),1) + (1:round(Fs*diff(timeRange)));

    % Set up multitaper
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = freqRange;
    winSize         = 0.25;
    winStep         = 0.025; % 4Hz resolution
    movingwin       = [winSize winStep];
    params.trialave = 1;

    e = load(fullfile(folderSegment,'LFP','elec1.mat'));
    goodTrials = setdiff(1:size(e.analogData,1),badTrials);
    numTrials = length(goodTrials);

    if numTrials>0
        for i=1:numElectrodes
            e = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeList(i)) '.mat']));
            [tfPower(i,:,:),timeValsTF0,freqVals] = mtspecgramc(e.analogData(goodTrials,goodTimePos)',movingwin,params); % for all the trials
            timeValsTF = timeValsTF0 + timeRange(1);
        end
    else
        tfPower = [];
        freqVals = [];
        timeValsTF = [];
    end
else
    disp([timingFile ' does not exist.']);
    tfPower = [];
    freqVals = [];
    timeValsTF = [];
    numTrials = 0;
end
end