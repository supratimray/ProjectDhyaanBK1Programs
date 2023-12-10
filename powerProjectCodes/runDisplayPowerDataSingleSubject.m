clear; close all

[allSubjectNames,expDateList] = getDemographicDetails('BK1');
[goodSubjectList, meditatorList, controlList] = getGoodSubjectsBK1;
folderSourceString = 'N:\Projects\ProjectDhyaan\BK1';
saveFolderName = 'powerResultsSingleSubject';

saveFileFlag     = 1;
plotRawTFFlag    = 0;
sortByBadTrialFlag = 1;

badEyeCondition = 'wo'; % use 'wo' for without, 'ep' for eye position and 'em' for eye movement
badTrialVersion = 'v8';
badTrialNameStr = ['_' badEyeCondition '_' badTrialVersion];

%%%%%%%%%%%%%%%%%%%%%% Bad trial rejection criteria %%%%%%%%%%%%%%%%%%%%%%%
badElectrodeRejectionFlag = 2; 
% 1: Don't reject badElectrodes,
% 2: reject badElectrodes for that protocol,
% 3: Reject badElectrodes common across all the protocols

useTheseIndices = 1:length(goodSubjectList);

for i=1:length(useTheseIndices)
    fh=figure(1); clf(fh);
    fh.WindowState = 'maximized';
    subjectName = goodSubjectList{useTheseIndices(i)};
    disp(['Analyzing for the subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};
    displayPowerDataSingleSubject(subjectName,expDate,folderSourceString,badTrialNameStr,badElectrodeRejectionFlag,plotRawTFFlag,sortByBadTrialFlag);

    if saveFileFlag
        makeDirectory(saveFolderName);
        fileNameTif = fullfile(saveFolderName,[subjectName badTrialNameStr '_badElecChoice' num2str(badElectrodeRejectionFlag) '_raw' num2str(plotRawTFFlag) '_sort' num2str(sortByBadTrialFlag) '.tif']);
        print(fh,fileNameTif,'-dtiff','-r300');
    end
end