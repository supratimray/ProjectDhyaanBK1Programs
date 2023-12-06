clear; close all

[allSubjectNames,expDateList] = getDemographicDetails('BK1');
[goodSubjectList, meditatorList, controlList] = getGoodSubjectsBK1;
folderSourceString      ='D:\Projects\ProjectDhyaan\BK1';

saveFileFlag     = 1;
closeFig         = 1;
badElectrodeList = [];
plotRawTFFlag    = [];

refScheme        = 1; % 1 for unipolar; 2 for bipolar
timeAvgFlag      = 1;
tfFlag           = 1;

allIndices = 1:length(goodSubjectList);
problamaticIndices   = find(strcmp(goodSubjectList,'053DR')); % this index need to be extracted separately
segmentTheseIndices  = setdiff(allIndices,problamaticIndices);
badTrialNameStr      ='_wo_v8';
badElectrodeRejectionFlag = 4; % 1: Don't reject badElectrodes,
% 2: reject badElectrodes for that protocol,
% 3: Reject badElectrodes common across all the protocols
% 4: Reject only the declared badElectrodes and reject protocol specific badElectrode

for i=1:length(segmentTheseIndices)
    fh=figure(1);
    fh.WindowState = 'maximized';
    subjectName = goodSubjectList{segmentTheseIndices(i)};
    disp(['Analyzing for the subject ' subjectName]);
    expDate = expDateList{strcmp(subjectName,allSubjectNames)};
    displayMeditationData(subjectName,expDate,folderSourceString,badTrialNameStr,badElectrodeRejectionFlag,plotRawTFFlag,refScheme,timeAvgFlag,tfFlag);

    if saveFileFlag
        fileName = [subjectName '_DisplayResultsV3' '.fig'];
        fileNameTif = [subjectName '_DisplayResultsV3'];
        saveas(fh,fileName);
        print(fh,fileNameTif,'-dtiff','-r300');
        if closeFig
            close(fh);
        end
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