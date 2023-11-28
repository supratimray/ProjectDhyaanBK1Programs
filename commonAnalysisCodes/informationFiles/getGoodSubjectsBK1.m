function [allSubjectList, meditatorList, controlList] = getGoodSubjectsBK1

declaredBadSubjects = [{'004P'} {'081SN'}];
tmp = load('BK1AllSubjectList.mat');
allSubjectList = setdiff(tmp.allSubjectList,declaredBadSubjects,'stable');
meditatorList = setdiff(tmp.meditatorList,declaredBadSubjects,'stable');
controlList = setdiff(tmp.controlList,declaredBadSubjects,'stable');
end