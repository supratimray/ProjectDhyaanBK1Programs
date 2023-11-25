
% input data format: 
%  n x 2 matrix. A 'violin' is plotted for each column 

% load the data
close all
clear
load("badElecNumberAllSub.mat");

dataArray{1,1} = badElecPercentageIndMeditators;
dataArray{1,2} = badElecPercentageIndConrols;

colorArray = {'r','g'};
showData = 1;
plotQuartiles = 1;
showSignificance = 1;
pairedDataFlag = 0;

figHandles=displayViolinPlot(dataArray,colorArray,showData,plotQuartiles,showSignificance,pairedDataFlag);