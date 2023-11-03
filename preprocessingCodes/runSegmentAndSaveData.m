% runSegmentAndSaveData
pairedSubjectNameList = getMatchedSubjectsBK1;
[subjectNameList,expDateList] = getDemographicDetails('BK1');

goodSubjectList = pairedSubjectNameList(:);

folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\ProjectDhyaan\BK1';

segmentTheseIndices = 13:length(goodSubjectList);
gridType = 'EEG'; 
capType = 'actiCap64_2019';
impedanceTag = '_Impedance_Start'; 
displayFlag = 0;
FsEye = 1000;

for i=1:length(segmentTheseIndices)
    subjectName = goodSubjectList{segmentTheseIndices(i)};
    disp(['Extracting Subject ' subjectName]);
    
    expDate = expDateList{strcmp(subjectName,subjectNameList)};
    segmentAndSaveData(subjectName,expDate,folderSourceString,FsEye); % Segment data
    getImpedanceDataEEG(subjectName,expDate,folderSourceString,gridType,impedanceTag,displayFlag,capType); % Get Impedance
end

% Errors found in
% i=5; 095KM
% i=12; 053DR