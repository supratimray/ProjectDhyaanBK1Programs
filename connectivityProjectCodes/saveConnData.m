
function saveConnData(subjectName,protocolNameList,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange)

if ~exist('stRange','var');         stRange = [0.25 1.25];              end

numProtocols = length(protocolNameList);
for i=1:numProtocols
    protocolName = protocolNameList{i};
    saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange);
end
end

function saveConnDataSingleProtocol(subjectName,protocolName,badEyeCondition,badTrialVersion,ftDataFolder,connMethod,stRange)

savedDataFolder = 'savedData'; % local project directory
makeDirectory(savedDataFolder);
savedDataFolder = fullfile(savedDataFolder,subjectName);
makeDirectory(savedDataFolder);

% Analysis file - this is where the final analyzed data is saved
analysisDetailsFileConn = fullfile(savedDataFolder,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']);

% ftData - data in fieldtrip format is saved in ftData
ftDataFileName = fullfile(ftDataFolder,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '.mat']);
tmpData = load(ftDataFileName);
numGoodTrials = tmpData.numGoodTrials;

if numGoodTrials>0    
    data = tmpData.data;

    cfg        = [];
    cfg.toilim = [(-diff(stRange) + 1/data.fsample) 0]; % just to have data of same length in Pre, Post
    dataPre   = ft_redefinetrial(cfg,data);
    
    cfg        = [];
    cfg.toilim = [(stRange(1) + 1/data.fsample) stRange(2)];
    dataPost   = ft_redefinetrial(cfg,data);

    [connPre,freqPre] = getConnIndividualSubject(dataPre,connMethod);
    [connPost,freqPost] = getConnIndividualSubject(dataPost,connMethod);
else
    connPre=[]; connPost=[]; freqPre=[]; freqPost=[];
end
save(analysisDetailsFileConn,'connPre','connPost','freqPre','freqPost','numGoodTrials');

end

function [conn,freqVals] = getConnIndividualSubject(data,connMethod)

temp=getConnThisCondition(data,connMethod);
conn = abs(temp.([connMethod 'spctrm']));
badElecs = data.badElecs;

conn(badElecs,:,:) = NaN;
conn(:,badElecs,:) = NaN;

freqVals = temp.freq;
end
function conn_stat = getConnThisCondition(data,method)
cfg              = [];
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.tapsmofrq    = 1;
if(strcmp(method,'coh'))
    cfg.output = 'powandcsd';
else
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
end
cfg.foilim       = [0 200];
freqPost = ft_freqanalysis(cfg, data);

if(strcmp(method,'coh'))
    TfreqPost = ft_checkdata(freqPost,'cmbrepresentation', 'full');
else
    TfreqPost = freqPost;
end

cfg=[];
cfg.method = method;
conn_stat = ft_connectivityanalysis(cfg, TfreqPost);
end