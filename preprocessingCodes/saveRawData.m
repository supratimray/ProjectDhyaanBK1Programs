% saveRawData
% The raw data obtained during the experiment is copied to an appropriate folder.

tmp = load('BK1AllSubjectList.mat');
allSubjectList = tmp.allSubjectList;

[subjectNameList,expDateList] = getDemographicDetails('BK1');

folderSourceString = 'N:\Projects\MeditationProjects\MeditationProject2\data\rawData';
folderDestinationString = 'N:\Projects\ProjectDhyaan\BK1\data\rawData';

for i=1:length(allSubjectList)
    subjectName = allSubjectList{i};
    expDate = expDateList{strcmp(subjectName,subjectNameList)};
    folderIn = fullfile(folderSourceString,[subjectName expDate]);
    folderOut = fullfile(folderDestinationString,[subjectName expDate]);

    if exist(folderOut,'file')
        disp([folderOut ' exists']);
    else
        disp(['Copying ' folderIn ' to ' folderOut]);
        copyfile(folderIn,folderOut,'f');
    end
end