% saveRawData
% The raw data obtained during the experiment is copied to an appropriate folder.

pairedSubjectNameList = getMatchedSubjectsBK1;
[subjectNameList,expDateList] = getDemographicDetails('BK1');

folderSourceString = 'N:\Projects\MeditationProjects\MeditationProject2\data\rawData';
folderDestinationString = 'N:\Projects\ProjectDhyaan\BK1\data\rawData';

for i=1:length(pairedSubjectNameList)
    for j=1:2
        subjectName = pairedSubjectNameList{i,j};
        expDate = expDateList{strcmp(subjectName,subjectNameList)};
        folderIn = fullfile(folderSourceString,[subjectName expDate]);
        folderOut = fullfile(folderDestinationString,[subjectName expDate]);

        disp(['Copying ' folderIn ' to ' folderOut]);
        copyfile(folderIn,folderOut,'f');
    end
end