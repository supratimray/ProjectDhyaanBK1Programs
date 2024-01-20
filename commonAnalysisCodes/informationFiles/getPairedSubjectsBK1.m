% 1st column: Meditator, 2nd Column: Control

function pairedSubjectNameList = getPairedSubjectsBK1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Males
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Male 20-29 
pairedSubjectNameList{1,1} = '019CKa'; pairedSubjectNameList{1,2} = '022SSP'; %  A-23, C-22 % Other options: 027SM not chosen because recording date is closer for 022SSP
pairedSubjectNameList{2,1} = '096MS';  pairedSubjectNameList{2,2} = '026HM';  %  A-26, C-26 % Other options: 021PB not chosen because of age (25 years)
pairedSubjectNameList{3,1} = '040VS';  pairedSubjectNameList{3,2} = '100UK';  %  A-27, C-28 % Other options: 026HM used earlier. 093AK better education matched with 012GK below
pairedSubjectNameList{4,1} = '012GK';  pairedSubjectNameList{4,2} = '093AK';  %  A-28, C-28 % Other options: 075AD has higher age difference. 100UK already used.
pairedSubjectNameList{5,1} = '095KM';  pairedSubjectNameList{5,2} = '075AD';  %  A-30, C-29 % Other options: 092KB is better matched education wise, but is the only match for meditator 051RA below.

% Male 30-39
% pairedSubjectNameList{6,1} = '051RA';  pairedSubjectNameList{6,2}  = '092KB'; % A-31, C-31. Removed now since 092KB has more than 24 bad electrodes and hence is decared as a bad subject
pairedSubjectNameList{6,1} = '090AV';  pairedSubjectNameList{6,2}  = '028HB'; % A-32, C-33 % Other options: 092KB used above
pairedSubjectNameList{7,1} = '054MP';  pairedSubjectNameList{7,2}  = '043AK'; % A-35, C-34 % Other options: 098GS is better age matched but is used for 015RK below
pairedSubjectNameList{8,1} = '015RK';  pairedSubjectNameList{8,2}  = '098GS'; % A-37, C-35 % Specal case: 098GS (35,M,13,0) used even though age difference is -2 years because otherwise there is no control
pairedSubjectNameList{9,1} = '038DK'; pairedSubjectNameList{9,2} = '071GK'; % A-38, C-39

% Male 40-44 
pairedSubjectNameList{10,1} = '041AG'; pairedSubjectNameList{10,2} = '077LK';  % A-41, C-40 % Other options: 003S has more difference in education level
pairedSubjectNameList{11,1} = '053DR'; pairedSubjectNameList{11,2} = '003S';   % A-43, C-43 % Other options: 063VK is better 'age-matched' but has too many bad electrodes 

pairedSubjectNameList{12,1} = '045SP';  pairedSubjectNameList{12,2} = '078BM'; % A-47, C-48
pairedSubjectNameList{13,1} = '025RK';  pairedSubjectNameList{13,2} = '070TB'; % A-49, C-49 % Other options: 078BM used earlier. 080RP has higher age differemce

% Male 50-59 
pairedSubjectNameList{14,1} = '035SS'; pairedSubjectNameList{14,2} = '080RP'; % A-50, C-50 % Other options: 048RU(51,M,NaN,0) and 070TB(49,M,24,0) have higher age difference
%pairedSubjectNameList{16,1} = '010AK'; pairedSubjectNameList{16,2} = '069MG'; % A-54, C-55. Removed now since 069MG has more than 24 bad electrodes and hence is decared as a bad subject
pairedSubjectNameList{15,1} = '044PN'; pairedSubjectNameList{15,2} = '101PB'; % A-54, C-52 % The only option is 069MG who is already taken. We again make an exception and take 101PB who is 52 years old.

% Male 60-65 
pairedSubjectNameList{16,1} = '046ME'; pairedSubjectNameList{16,2} = '085BM'; % A-62, C-61 % Other option: 068RV whos education level is not known

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Females
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Female 20-29 
pairedSubjectNameList{17,1} = '031BK'; pairedSubjectNameList{17,2} = '073SK'; % A-24, C-26; % '081SN' is a better match age wise (24 years) but has too many bad electrodes and hence not used
pairedSubjectNameList{18,1} = '056PR'; pairedSubjectNameList{18,2} = '086AB'; % A-27, C-28: Other option: 073SK is used above.

% Female 30-39
pairedSubjectNameList{19,1} = '052PR'; pairedSubjectNameList{19,2} = '082MS'; % A-31,  C-31 % Other option: 087KR has smaller difference in MC phase but has too many bad eye trials for G2
pairedSubjectNameList{20,1} = '059MS'; pairedSubjectNameList{20,2} = '102AS'; % A-34,  C-35 % Other option: 064PK has larger difference in Education level
pairedSubjectNameList{21,1} = '013AR'; pairedSubjectNameList{21,2} = '064PK'; % A-35,  C-35 % Other option: 102AS used above
pairedSubjectNameList{22,1} = '074KS'; pairedSubjectNameList{22,2} = '084AK'; % A-36,  C-37 % Other options: 064PK and 102AS used above

% Female 40-49  
pairedSubjectNameList{23,1}  = '006SR';  pairedSubjectNameList{23,2}  = '062MT'; % A-41, C-42
pairedSubjectNameList{24,1}  = '094SR';  pairedSubjectNameList{24,2}  = '049KK'; % A-45, C-46
pairedSubjectNameList{25,1}  = '017KG';  pairedSubjectNameList{25,2}  = '072DK'; % A-49, C-49

% Female 50-59 
pairedSubjectNameList{26,1} = '042VA'; pairedSubjectNameList{26,2} = '097SV'; % A-50, C-51 % Other options: 072DK and 099SP - both have higher difference in education level
pairedSubjectNameList{27,1} = '060GV'; pairedSubjectNameList{27,2} = '099SP'; % A-51, C-51 
pairedSubjectNameList{28,1} = '050UR'; pairedSubjectNameList{28,2} = '083SP'; % A-53, C-53
pairedSubjectNameList{29,1} = '030SH'; pairedSubjectNameList{29,2} = '076BH'; % A-56, C-56

% Female 60-65 
pairedSubjectNameList{30,1} = '089AB'; pairedSubjectNameList{30,2} = '079SG'; % A-62, C-62
end
