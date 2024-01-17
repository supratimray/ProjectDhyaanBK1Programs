% display TF data for Matched Subjects (A and C)

clear
close all

protocolIndex = 3;
numMatchedMales = 18;
numMatchedFemales = 14;
protocolNameList = [{'G1'} {'G2'} {'M2'}];

diffTf = 1;
fontsize = 10;
badElecRejectionFlag = 1;
baselineRange = [-1 0];

slowGammaRange = [22 34];
fastGammaRange = [35 65];

% getting the demographic information
pairedSubjectNameList = getPairedSubjectsBK1;
numPairs  = length(pairedSubjectNameList);
numGroups = size(pairedSubjectNameList,2);
[subjectNameList,expDateList,labelList,ageList,genderList,educationList,mcList] = getDemographicDetails('BK1');

% get the agesList
ageListPaired = zeros(numPairs,2);
for j=1:numGroups
    for i=1:numPairs
        subIndex = find(strcmp(subjectNameList,pairedSubjectNameList(i,j)));
        ageListPaired(i,j) = ageList(subIndex);
    end
end

% plot for a single subject:
hTF1 = getPlotHandles(numMatchedMales,numGroups,[0.08 0.06 0.4 0.84],0.005,0.005);
hTF2 = getPlotHandles(numMatchedMales,numGroups,[0.54  0.06 0.4 0.84],0.005,0.005);
hTF = [hTF1;hTF2];

for i=1:numPairs
    for j=1:numGroups
        subjectName = pairedSubjectNameList{i,j};
        disp(subjectName);
        fileName = fullfile(pwd,'savedData',[subjectName '_ep_v8_TF.mat']);
        load(fileName);

        if badElecRejectionFlag
            badElecToReject = badElectrodes{protocolIndex};
            goodElecsInd    = not(ismember(electrodeList,badElecToReject));
            if sum(goodElecsInd)>1
                meanTfPower = squeeze(mean(tfPower{protocolIndex}(goodElecsInd,:,:),1));
            else
                disp('Not enough good electrodes to plot')
                continue
            end
        else
            meanTfPower = squeeze(mean(tfPower{protocolIndex}(:,:,:),1));
        end

        logP = log10(meanTfPower);
        baselinePower = mean(logP(timeValsTF>=baselineRange(1) & timeValsTF<=baselineRange(2),:));
        if diffTf
            pcolor(hTF(i,j),timeValsTF,freqValsTF,10*(logP'- repmat(baselinePower',1,length(timeValsTF))));
        else
            pcolor(hTF(i,j),timeValsTF,freqValsTF,logP');
        end

        % Add labels
        axes(hTF(i,j));
        shading(hTF(i,j),'interp'); colormap jet;
        yline(slowGammaRange(1),'k--');
        yline(slowGammaRange(2),'k--');
        text(-0.2,50,[num2str(ageListPaired(i,j)) 'Yrs'],'Parent',hTF(i,j));
        xlim(hTF(i,j),[-0.25 max(timeValsTF)]);
        set(hTF(i,j), 'TickDir', 'out');

        if ismember(i,[numMatchedMales,numPairs]) && j==1
            xlabel(hTF(i,j),'Time(s)');
            ylabel(hTF(i,j),'Frequency (Hz)');
            set(hTF(i,j), 'fontsize',fontsize);
            set(hTF(i,j), 'FontWeight','bold');
        else
            set(hTF(i,j),'Yticklabel',[]);
            set(hTF(i,j),'Xticklabel',[]);
        end

        if j==1 && ismember(i,[1,numMatchedMales+1])
            if i==1
                title(hTF(i,j),'Meditators(M)');
            else
                title(hTF(i,j),'Meditators(F)');
            end
        elseif j==2 && ismember(i,[1,numMatchedMales+1])
            if i==1
                title(hTF(i,j),'Controls(M)');
            else
                title(hTF(i,j),'Controls(F)');
            end
        end

        clim(hTF(i,j),[-5 5]);
        ylim(hTF(i,j),[0 70]);
    end
end

% set the colorbar
set(hTF(i,j), 'TickDir', 'out');
hc = colorbar('Position', [0.05 0.2 0.02 0.3]);

hc.FontSize         = 10;
hc.Label.FontSize   = 10;
hc.Label.FontWeight = 'bold';
hc.Label.String = ['\Delta Power' '(dB)'];

% hide remaining plots:
numHidePlots = numMatchedMales-numMatchedFemales;
set(hTF(numPairs+1:numPairs+numHidePlots,1:2), 'visible', 'off');
sgtitle(['TF plot for Med vs Control (' protocolNameList{protocolIndex} ')']);
