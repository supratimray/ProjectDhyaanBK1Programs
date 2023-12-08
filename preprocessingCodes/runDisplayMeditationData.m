clear; close all

[allSubjectNames,expDateList] = getDemographicDetails('BK1');
[goodSubjectList, meditatorList, controlList] = getGoodSubjectsBK1;
folderSourceString      ='N:\Projects\ProjectDhyaan\BK1';
saveFolderName = 'Results';

saveFileFlag     = 1;
plotRawTFFlag    = 0;

badEyeCondition = 'wo'; % use 'wo' for without, 'ep' for eye position and 'em' for eye movement
badTrialVersion = 'v8';
badTrialNameStr = ['_' badEyeCondition '_' badTrialVersion];

%%%%%%%%%%%%%%%%%%%%%% Bad trial rejection criteria %%%%%%%%%%%%%%%%%%%%%%%
badElectrodeRejectionFlag = 3; 
% 1: Don't reject badElectrodes,
% 2: reject badElectrodes for that protocol,
% 3: Reject badElectrodes common across all the protocols

useTheseIndices = 1;

for i=1:length(useTheseIndices)
    fh=figure(1); clf(fh);
    fh.WindowState = 'maximized';
    subjectName = goodSubjectList{useTheseIndices(i)};
    disp(['Analyzing for the subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};
    displayMeditationData(subjectName,expDate,folderSourceString,badTrialNameStr,badElectrodeRejectionFlag,plotRawTFFlag);

    if saveFileFlag
        makeDirectory(saveFolderName);
        fileNameTif = fullfile(saveFolderName,[subjectName badTrialNameStr '_badElecChoice' num2str(badElectrodeRejectionFlag) '_raw' num2str(plotRawTFFlag) '.tif']);
        print(fh,fileNameTif,'-dtiff','-r300');
    end
end

% Issues during the data display process:
% Fixed by putting extra conditions while plottng!
% i = 21, 059MS; EC2 segment
% i = 34; 024SK; G2 runs for 180  trials; 
% i = 45; 063VK; Group 4 
% i = 58; 080RP
% i = 66; 092KB
% i = 67; 093AK
% i = 69; 098GS;
% i = 70; 099SP; 