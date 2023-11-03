%%%%%%%%%%%%%%%%%%%%%% Get Subject Details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[subjectNameList,expDateList,labelList,ageList,genderList,educationList,mcList] = getDemographicDetails('BK1');

meditatorPos = strcmp(labelList,'A');
beginnerPos = strcmp(labelList,'B');
controlPos =  strcmp(labelList,'C');
malePos = strcmp(genderList,'M');
femalePos = strcmp(genderList,'F');

disp(['Beginners (discarded): ' num2str(length(find(beginnerPos)))]);
disp(['Advanced (meditator): ' num2str(length(find(meditatorPos))) '(M=' num2str(length(find(meditatorPos & malePos))) ',F=' num2str(length(find(meditatorPos & femalePos))) ')']);
disp(['Control: ' num2str(length(find(controlPos))) '(M=' num2str(length(find(controlPos & malePos))) ',F=' num2str(length(find(controlPos & femalePos))) ')']);

%%%%%%%%%%%%%%%%%%%%%% Find Matched subjects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1 - For each meditator, find controls of same gender and age within +- ageLim years
ageLim = 1;

% Create meditator list sorted by gender and age
maleMeditatorPos = find(meditatorPos & malePos);
maleMeditatorAge = ageList(maleMeditatorPos);
[~,order] = sort(maleMeditatorAge);
meditatorList = subjectNameList(maleMeditatorPos(order));

femaleMeditatorPos = find(meditatorPos & femalePos);
femaleMeditatorAge = ageList(femaleMeditatorPos);
[~,order] = sort(femaleMeditatorAge);
meditatorList = cat(1,meditatorList,subjectNameList(femaleMeditatorPos(order)));

for i=1:length(meditatorList)
    subjectName = meditatorList{i};
    pos = find(strcmp(subjectName,subjectNameList));
    age = ageList(pos); gender = genderList(pos);
    
    ageMatchPos = (ageList<= age+ageLim) & (ageList>= age-ageLim);
    genderMatchPos = strcmp(gender,genderList);
    allPos = controlPos & ageMatchPos & genderMatchPos;
    controls = subjectNameList(allPos);
    
    controlStr = '';
    for j=1:length(controls)
        cPos = find(strcmp(controls{j},subjectNameList));
        controlStr = cat(2,controlStr,[' ' controls{j} ' (' num2str(ageList(cPos)) ',' genderList{cPos} ',' num2str(educationList(cPos)) ',' num2str(mcList(cPos))  ')']);
    end
    disp([num2str(i) ': ' subjectName ' (' num2str(age) ',' gender{1} ',' num2str(educationList(pos)) ',' num2str(mcList(pos))  '). Control(s): ' controlStr]);
end
