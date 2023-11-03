
function [subjectNameList,expDateList,labelList,ageList,genderList,educationList,mcList] = getDemographicDetails(projectName)

d = load([projectName 'ProjectDetails.mat']);
demographicDetails = d.demographicDetails;
N = length(demographicDetails);

[~,nameRow] = find(strcmpi(demographicDetails,'SubjectName'));
subjectNameList = demographicDetails(2:N,nameRow);

[~,expDateRow] = find(strcmpi(demographicDetails,'ExpDate'));
expDateList = demographicDetails(2:N,expDateRow);

[~,labelRow] = find(strcmpi(demographicDetails,'Label'));
labelList = demographicDetails(2:N,labelRow);

[~,ageRow] = find(strcmpi(demographicDetails,'Age'));
ageList = cell2mat(demographicDetails(2:N,ageRow));

[~,genderRow] = find(strcmpi(demographicDetails,'Gender'));
genderList = demographicDetails(2:N,genderRow);

[~,educationRow] = find(strcmpi(demographicDetails,'Education in Years'));
educationList = cell2mat(demographicDetails(2:N,educationRow));

[~,mcRow] = find(strcmpi(demographicDetails,'Menstrual Cycle'));
mcList = cell2mat(demographicDetails(2:N,mcRow));

end