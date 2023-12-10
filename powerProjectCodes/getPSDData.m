function [psdValsST,psdValsBL,freqVals,numTrials,badElectrodes] = getPSDData(subjectName,expDate,protocolNameList,folderSourceString,badEyeCondition,badTrialVersion,stRange)

if ~exist('protocolNameList','var')
    protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}]; 
end
if ~exist('folderSourceString','var');   folderSourceString=[];         end
if ~exist('badEyeCondition','var');      badEyeCondition = 'wo';        end
if ~exist('badTrialVersion','var');      badTrialVersion = 'v8';        end
if ~exist('stRange','var');              stRange = [0.25 1.25];         end
if isempty(folderSourceString)
    folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Protocols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numProtocols = length(protocolNameList);

psdValsST = cell(1,numProtocols);
psdValsBL = cell(1,numProtocols);
numTrials = zeros(1,numProtocols);

badElectrodes = cell(1,numProtocols);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    [badTrials,badElectrodes{i}] = getBadTrialsAndElectrodes(subjectName,expDate,protocolName,folderSourceString,badEyeCondition,badTrialVersion);
    [psdValsST{i},psdValsBL{i},freqVals,numTrials(i)] = getDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,stRange,badTrials);
end
end

function [psdValsST,psdValsBL,freqVals,numTrials] = getDataSingleProtocol(subjectName,expDate,protocolName,folderSourceString,stRange,badTrials)

gridType = 'EEG'; 
electrodeList = 1:64;
tapers = [1 1];
freqRange = [0 150];

numElectrodes = length(electrodeList);

folderSegment = fullfile(folderSourceString,'data','segmentedData',subjectName,gridType,expDate,protocolName,'segmentedData');
timingFile = fullfile(folderSegment,'LFP','lfpInfo.mat');

if exist(timingFile,'file')

    t = load(timingFile);
    timeVals = t.timeVals;
    Fs = round(1/(timeVals(2)-timeVals(1)));
    goodTimePosST = find(timeVals>=stRange(1),1) + (1:round(Fs*diff(stRange)));

    blRange = [-diff(stRange) 0];
    goodTimePosBL = find(timeVals>=blRange(1),1) + (1:round(Fs*diff(blRange)));

    % Set up multitaper
    params.tapers   = tapers;
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = freqRange;
    params.trialave = 1;

    e = load(fullfile(folderSegment,'LFP','elec1.mat'));
    goodTrials = setdiff(1:size(e.analogData,1),badTrials);
    numTrials = length(goodTrials);

    if numTrials>0
        for i=1:numElectrodes
            e = load(fullfile(folderSegment,'LFP',['elec' num2str(electrodeList(i)) '.mat']));
            [psdValsST(i,:),freqVals] = mtspectrumc(e.analogData(goodTrials,goodTimePosST)',params); %#ok<*AGROW>
            psdValsBL(i,:) = mtspectrumc(e.analogData(goodTrials,goodTimePosBL)',params);
        end
    else
        psdValsST = [];
        psdValsBL = [];
        freqVals = [];
    end
else
    disp([timingFile ' does not exist.']);
    psdValsST = [];
    psdValsBL = [];
    freqVals = [];
    numTrials = 0;
end
end