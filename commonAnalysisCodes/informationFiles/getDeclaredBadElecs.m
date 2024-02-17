function [declaredBadElectrodes] = getDeclaredBadElecs(threshold)
% Gives the declared bad electrodelist according to the current threshold
% We get the this list after setting the desired thresold in the runDisplayBadElecsAcrossSubjectsBK1 script

if ~exist('threshold','var');   threshold=0.35;       end

if threshold==0.35
    declaredBadElectrodes = [1 5 28 38 39 55 57 58]; % threshold 0.35
end
end

