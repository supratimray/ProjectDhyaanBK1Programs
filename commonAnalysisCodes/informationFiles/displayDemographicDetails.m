function displayDemographicDetails(projectName)

[subjectNameList,~,labelList,ageList,genderList,educationList,mcList] = getDemographicDetails(projectName);

% Get Entries for each usable subject

genderIDs = [{'M'} {'F'}];
labelIDs = [{'A'} {'B'} {'C'}];

% Create table
numColumns = 4;
ageCutoff = 74;
numRows = ageCutoff-18; % Age 20-65 + 1 row for title

hPlots = getPlotHandles(1,2,[0.05 0.05 0.9 0.9]);

for g=1:2
    tableContents = cell(numRows,numColumns);
    tableContents{1,1} = 'Age';
    posToUse = find(strcmp(genderList,genderIDs{g}));
    labelIDsToUse = labelList(posToUse);
    for i=1:3
        tableContents{1,i+1} = [labelIDs{i} ' (' num2str(length(find(strcmp(labelIDsToUse,labelIDs{i})))) ')'];
    end

    for i=1:length(posToUse)
        pos = posToUse(i);
        age = ageList(pos);
        groupNum = find(strcmp(labelIDs,labelList(pos)));

        rowNum = age-18;
        colNum = groupNum+1;
        tableContents{rowNum,1} = num2str(age);
        if isempty(tableContents{rowNum,colNum})
            if g==1
                tableContents{rowNum,colNum} = [tableContents{rowNum,colNum} [subjectNameList{pos} '(' num2str(educationList(pos)) ')']];
            else
                tableContents{rowNum,colNum} = [tableContents{rowNum,colNum} [subjectNameList{pos} '(' num2str(educationList(pos)) ',' num2str(mcList(pos)) ')']];
            end
        else
            if g==1
                tableContents{rowNum,colNum} = [tableContents{rowNum,colNum} '/' [subjectNameList{pos} '(' num2str(educationList(pos)) ')']];
            else
                tableContents{rowNum,colNum} = [tableContents{rowNum,colNum} '/' [subjectNameList{pos} '(' num2str(educationList(pos)) ',' num2str(mcList(pos)) ')']];
            end
        end
    end

    makeTable(hPlots(g),tableContents);
    title(hPlots(g),[genderIDs{g} '(' num2str(length(posToUse)) ')']);
end
end