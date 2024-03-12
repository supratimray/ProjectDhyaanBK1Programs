% subjectNameLists - 1x2 cell array containing subjectNames for meditators and controls
% protocolName - one of EO1, EC1, G1, M1, G2, EC2, EO2, M2
% analysisChoice - 'st', 'bl' or 'combined'
% ElectrodeList - electrodes for which connectivity should be measured and averaged
% connMethod - connectivity method
% Option added to return connectivity and topoplot data. Also to simply return these without displaying here.

function [connDataToReturn,topoplotDataToReturn,freqVals,binnedCenters] = displayConnDataAllSubjects(subjectNameLists,protocolName,analysisChoice,electrodeList,connMethod,badEyeCondition,badTrialVersion,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,displayDataFlag)

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('electrodeList','var');         electrodeList = [];           end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [3 30];
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

numFreqRanges = length(freqRangeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [0.8 0 0.8];      % Purple 
displaySettings.colorNames(2,:) = [0.25 0.41 0.88]; % Cyan
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

cLimsTopo = axisRangeList{2};
binRange = [-0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%% Get topoplot info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capType = 'actiCap64_UOL';
x = load([capType '.mat']); 
montageChanlocs = x.chanlocs;

saveFolderName = 'savedData';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    if isempty(hAllPlots)
        hTopo = getPlotHandles(numFreqRanges,2,[0.05 0.05 0.4 0.9],0.025,0.025,1);
        hConn = getPlotHandles(numFreqRanges,1,[0.5 0.05 0.2 0.9],0.025,0.025,0);
        hBar  = getPlotHandles(numFreqRanges,1,[0.75 0.05 0.15 0.9],0.025,0.025,0);
    else
        hTopo = hAllPlots.hTopo;
        hConn = hAllPlots.hConn;
        hBar  = hAllPlots.hBar;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[connData,freqVals,connDataBinwise,binnedCenters] = getConnDataAllSubjects(subjectNameLists,electrodeList,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,montageChanlocs);

%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
lineNoiseRange = [48 52];
badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
for i = 1:numFreqRanges
    freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<freqRangeList{i}(2))),badFreqPos);
end

goodBinPos = intersect(find(binnedCenters>=binRange(1)),find(binnedCenters<=binRange(2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
topoplotDataToReturn = cell(2,numFreqRanges);
connDataToReturn = cell(2,numFreqRanges);
dataBar = cell(2,numFreqRanges);

for i=1:2
    for j=1:numFreqRanges

        % Topoplots
        x = squeeze(mean(connData{i}(:,:,freqPosList{j}),3));
        if useMedianFlag
            data = squeeze(median(x,1,'omitnan'));
        else
            data = squeeze(mean(x,1,'omitnan'));
        end
        
        topoplotDataToReturn{i,j} = data;
        if displayDataFlag
            axes(hTopo(j,i)); %#ok<*LAXES>
            topoplot(data,montageChanlocs,'electrodes','on','maplimits',cLimsTopo,'plotrad',0.6,'headrad',0.6); colorbar;
        end

        % Conn versus distance
        x = squeeze(mean(connDataBinwise{i}(:,:,freqPosList{j}),3));
        connDataToReturn{i,j} = x;
        if useMedianFlag
            data = squeeze(median(x,1,'omitnan'));
        else
            data = squeeze(mean(x,1,'omitnan'));
        end
        plot(hConn(j),binnedCenters,data,'color',displaySettings.colorNames(i,:));
        hold(hConn(j),'on');

        % barplot
        dataBar{i,j} = squeeze(mean(x(:,goodBinPos),2,'omitnan'));
    end
end

% show bar plots
if displayDataFlag
    for i=1:numFreqRanges

        % display violin plots for power
        displaySettings.plotAxes = hBar(i);
        if ~useMedianFlag
            displaySettings.parametricTest = 1;
            displaySettings.medianFlag = 0;
        else
            displaySettings.parametricTest = 0;
            displaySettings.medianFlag = 1;
        end

        displaySettings.showYTicks=1;
        displaySettings.showXTicks=1;
        displaySettings.yPositionLine=0.5;
        displaySettings.setYLim = [0 1];

        displayViolinPlot(dataBar(:,i)',[{displaySettings.colorNames(1,:)} {displaySettings.colorNames(2,:)}],1,1,1,pairedDataFlag,displaySettings);
        ylim(hBar(i),[0 1]);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [connData,freqVals,connDataBinwise,binnedCenters] = getConnDataAllSubjects(subjectNameLists,electrodeList,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,montageChanlocs)

%%%%% Discretize connectivty into bins depending on distance from seed %%%%
binWidth = 0.25;
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binnedCenters = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(montageChanlocs);

numElectrodes = length(electrodeList);
binnedIndicesAllElectrodes = cell(numElectrodes,nbins);

for e=1:numElectrodes
    dist = sqrt((angl_dist(loc.azi(electrodeList(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(electrodeList(e)),loc.ele,'e')).^2);
    fitx = cos((dist/180)*pi);
    binned_fitx = discretize(fitx,binEdges);
    
    for b = 1:nbins % number of bins
        binnedIndicesAllElectrodes{e,b} = find(binned_fitx == b);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
badSubjectList = cell(1,2);
connDataTMP = cell(1,2);
connDataBinwiseTMP = cell(1,2);

for i=1:2
    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod]));
        numGoodTrials = tmpData.numGoodTrials;

        if numGoodTrials<cutoffList(2)
            badSubjectList{i}(j) = 1;
        else
            if strcmp(analysisChoice,'bl')
                connDataTMP2 = tmpData.connPre(electrodeList,:,:);
                freqVals = tmpData.freqPre;
            elseif strcmp(analysisChoice,'st')
                connDataTMP2 = tmpData.connPost(electrodeList,:,:);
                freqVals = tmpData.freqPost;
            else
                connDataTMP2 = (tmpData.connPre(electrodeList,:,:) + tmpData.connPost(electrodeList,:,:))/2;
                freqVals = tmpData.freqPost;
            end

            numGoodElectodes = trace(~isnan(squeeze(connDataTMP2(:,electrodeList,1))));

            if numGoodElectodes >= cutoffList(1)
                badSubjectList{i}(j) = 0;
                connDataTMP{i}{j} = squeeze(mean(connDataTMP2,1,'omitnan'));
                
                connDataBinwiseTMP2 = zeros(numElectrodes,nbins,length(freqVals));
                for e=1:numElectrodes
                    for b=1:nbins
                         connDataBinwiseTMP2(e,b,:) = squeeze(mean(connDataTMP2(e,binnedIndicesAllElectrodes{e,b},:),2,'omitnan'));
                    end
                end
                connDataBinwiseTMP{i}{j} = squeeze(mean(connDataBinwiseTMP2,1,'omitnan'));
            else
                badSubjectList{i}(j) = 1;
            end
        end
    end
end

% Remove bad subjects
connData = cell(1,2);
connDataBinwise = cell(1,2);

for i=1:2
    if pairedDataFlag
        badSubjectPos = find(sum(cell2mat(badSubjectList')));
    else
        badSubjectPos = find(badSubjectList{i});
    end
    x1 = connDataTMP{i};
    x1(badSubjectPos)=[];
    x2 = connDataBinwiseTMP{i};
    x2(badSubjectPos)=[];

    numSubjects = length(x1);
    for j=1:numSubjects
        connData{i}(j,:,:) = x1{j};
        connDataBinwise{i}(j,:,:) = x2{j};
    end
end
end

%%%%%%%% Codes written by Santosh for TLSA connectivity project %%%%%%%%%%%
function loc = getElecLocAngles(chanlocs)
azi = zeros(1,length(chanlocs)); ele = zeros(1,length(chanlocs));
for e = 1:length(chanlocs)
    azi(e) = chanlocs(e).sph_theta;
    ele(e) = chanlocs(e).sph_phi;
end
loc.azi = azi;
loc.ele = ele;
end
function out_theta = angl_dist(in_theta_ref,in_theta,val)
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    if(in_theta_ref > 90)
    in_theta(14) = 90;
    elseif(in_theta_ref < -90)
        in_theta(14) = -90;
    end
end
in_theta = abs(in_theta-in_theta_ref);
out_theta = zeros(1,length(in_theta));
for i=1:length(in_theta)
    if(in_theta(i) > 180)
        out_theta(i) = 360 - in_theta(i); % to get shortest angular distance
    else
        out_theta(i) = in_theta(i);
    end
end
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    out_theta(14) = 0;
end
end