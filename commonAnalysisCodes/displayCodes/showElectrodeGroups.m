function montageChanlocs = showElectrodeGroups(hPlots,capType,electrodeGroupList,groupNameList, displaySettings)

%%%%%%%%%%%%%%%%%%%%%% check the displaySettings %%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('displaySettings','var');                displaySettings=struct();              end
if ~isfield(displaySettings,'fontSizeEG');         displaySettings.fontSizeEG = 10;       end
if ~isfield(displaySettings,'headRad');            displaySettings.headRad = 0.6;         end
if ~isfield(displaySettings,'colormap');           displaySettings.colormap = 'jet';      end  

%%%%%%%%%%%%%%%%%%%%%% Compare with Montage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = load([capType 'Labels.mat']); montageLabels = x.montageLabels(:,2);
x = load([capType '.mat']);       montageChanlocs = x.chanlocs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Topoplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(hPlots(1));
electrodeSize = 5;
numElectrodes = length(montageLabels);

for i=1:numElectrodes
    montageChanlocs(i).labels = ' ';
end

numElectrodeGroups = length(electrodeGroupList);

if strcmp(displaySettings.colormap,'lines')
    electrodeGroupColorList = lines(numElectrodeGroups);
else
   electrodeGroupColorList = jet(numElectrodeGroups); 
end

for i=1:numElectrodeGroups
    topoplot(zeros(1,numElectrodes),montageChanlocs,'maplimits',[-1 1],'electrodes','on','style','map','emarker2',{electrodeGroupList{i},'o',electrodeGroupColorList(i,:),electrodeSize},'plotrad',0.6,'headrad',displaySettings.headRad);
end
topoplot([],montageChanlocs,'electrodes','labels','style','blank','plotrad',0.6,'headrad',displaySettings.headRad);

axes(hPlots(2))
set(hPlots(2),'visible','off');
for i=1:numElectrodeGroups
    text(0.05,0.9-0.15*(i-1),[groupNameList{i} ' (N=' num2str(length(electrodeGroupList{i})) ')'],'color',electrodeGroupColorList(i,:),'unit','normalized','Fontsize', displaySettings.fontSizeEG);
end
end