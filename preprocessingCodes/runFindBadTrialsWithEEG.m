% runFndBadTrialsWithEEG for the project BK1

% folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\ProjectDhyaan\BK1';
folderSourceString = 'D:\Projects\ProjectDhyaan\BK1';
[subjectNameList,expDateList] = getDemographicDetails('BK1');

% get the paired subjectList:
% pairedSubjectNameList = getMatchedSubjectsBK1;
% goodSubjectList = pairedSubjectNameList(:);

% get all the subjectList:
fileToLoadSourceString = 'D:\Projects\ProjectDhyaan\BK1\ProjectDhyaanBK1Programs\commonAnalysisCodes\informationFiles';
fileName = 'BK1AllSubjectList.mat';
load(fullfile(fileToLoadSourceString,fileName),'allSubjectList');

allIndices = 1:length(allSubjectList);
problamaticIndices  = 12; % this index need to be run separately with minor chnage (se below)
segmentTheseIndices = setdiff(allIndices,problamaticIndices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gridType = 'EEG';
nonEEGElectrodes = 65:70;
impedanceTag = '_Impedance_Start';
capType = 'actiCap64_UOL';
saveDataFlag = 1;

badTrialNameStrS = {'_wo_v8'};
displayResultsFlag = 1;
electrodeGroup = ''; % by default, highPriorityElectrodes are used
checkPeriod = [-1 1.25];
checkBaselinePeriod = checkPeriod;
badEEGElectrodes = [];

useEyeData =       [   0       0       0     0      0       0      0       0    ];
protocolNameList = [{'EO1'} {'EC1'} {'G1'} {'M1'} {'G2'} {'EO2'} {'EC2'} {'M2'}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b = 1:length(badTrialNameStrS)
    badTrialNameStr=badTrialNameStrS{b};
    for i=segmentTheseIndices
        subjectName = allSubjectList{i};
        disp(['Extracting Subject ' subjectName]);
        expDate = expDateList{strcmp(subjectName,subjectNameList)};
        for j=1:length(protocolNameList)
            protocolName = protocolNameList{j};
            useEyeDataFlag = useEyeData(j);
            findBadTrialsWithEEG(subjectName,expDate,protocolName,folderSourceString,gridType,badEEGElectrodes,...
                nonEEGElectrodes,impedanceTag,capType,saveDataFlag,badTrialNameStr,displayResultsFlag,electrodeGroup,checkPeriod,checkBaselinePeriod,useEyeDataFlag);
        end
    end
end

%-----------------------------------------------------------------------------------------------------
% Issues during finding the badtrials:
% problamaticIndices  = [10,12];
% goodSegmentIndices  = setdiff(segmentTheseIndices,problamaticIndices);
%------------------------------------------------------------------------------------------------------
% i=10; 048RU
% Impedance measurement was done only for the first 32 channels.
% We have made impedance values of the rest 32 channelss as zero during segementation
% i=12; 053DR
% EO1 missing; Bad trial code was run for the rest of the segments.
%-----------------------------------------------------------------------------------------------------