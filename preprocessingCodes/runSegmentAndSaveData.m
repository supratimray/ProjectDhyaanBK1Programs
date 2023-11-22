% runSegmentAndSaveData for the BK1 project

gridType = 'EEG';
capType = 'actiCap64_2019';
impedanceTag = '_Impedance_Start';
displayFlag = 0;
FsEye = 1000;
% folderSourceString = 'D:\OneDrive - Indian Institute of Science\Supratim\Projects\ProjectDhyaan\BK1';
folderSourceString = 'D:\Projects\ProjectDhyaan\BK1';
[subjectNameList,expDateList] = getDemographicDetails('BK1');

% get the paired subjectList
% pairedSubjectNameList = getMatchedSubjectsBK1;
% goodSubjectList = pairedSubjectNameList(:);

% get all the subjectList
fileToLoadSourceString = 'D:\Projects\ProjectDhyaan\BK1\ProjectDhyaanBK1Programs\commonAnalysisCodes\informationFiles';
fileName = 'BK1AllSubjectList.mat';
load(fullfile(fileToLoadSourceString,fileName),'allSubjectList');

allIndices = 1:length(allSubjectList);
problamaticIndices  = [5 12 75]; % this indices need to be extracted separately 
segmentTheseIndices  = setdiff(allIndices,problamaticIndices);

for i=segmentTheseIndices
    subjectName = allSubjectList{i};
    disp(['Extracting Subject ' subjectName]);
    %-------------------------------------------------------------------------------------------------------------------------------
    expDate = expDateList{strcmp(subjectName,subjectNameList)};
    segmentAndSaveData(subjectName,expDate,folderSourceString,FsEye); % Segment data
    getImpedanceDataEEG(subjectName,expDate,folderSourceString,gridType,impedanceTag,displayFlag,capType); % Get Impedance
end

%-----------------------------------------------------------------------------------------------------
% Issues during the data segment process:
% problamaticIndices  = [5 12 75];
% goodSegmentIndices  = setdiff(segmentTheseIndices,problamaticIndices);
% ----------------------------------------------------------------------------------------------------
% Marker Issues:-
% i=5;  095KM, M2 protocol, .vmrk file had junk markers; for extraction we have used the markers which matches with the .bhv2 file;
% i=75; 099SP, M2 protocol has less trial start markers compared to .bhv2 file. Saved data for only the trials which have markers matched between
%       .vmrk and .bhv2 file. 2 trials were lost at the beggining of the experiment
% ----------------------------------------------------------------------------------------------------
% File Missing:-
% i=12; 053DR, EO1 protocol could not be extracted as .bhv2 file is missing;
%-----------------------------------------------------------------------------------------------------
% Time issues:- (could be extracted with the updated version of the segmentation code)
% Fixed with additional condition while processing the eye data in the segmentation codes
% i=15; 035SS, EC2
% i=26; 006SR, EO1 protocol
% i=35; 089AB, EO1 and EC1 protocol
% i=50; 064PK, EO1 and EC1 protocol
%-----------------------------------------------------------------------------------------------------