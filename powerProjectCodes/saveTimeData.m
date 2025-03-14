function protcolLengthTime = saveTimeData(fileNameSegmented,fileNameSave,protocolNameMain,subjectName,appendFlag)

if ~exist('appendFlag','var'); appendFlag=0; end

protcolLengthTime = zeros(1,length(protocolNameMain));
for i=1:length(protocolNameMain)
    % read file
    protocolToProcess = protocolNameMain{i};
    fileNameRead = fullfile(fileNameSegmented,protocolToProcess,'extractedData','goodStimCodeNums.mat');
    % load file into a temp variable
    x=load(fileNameRead);
    times = x.goodStimTimes;
    protcolLengthTime(i) = round((times(end)-times(1))/60,2);
end
% append the data
if appendFlag
    save(fileNameSave, 'protcolLengthTime', '-append');
    fprintf('Time data saved successfully for subject %s\n', subjectName);
end
end