function [eyeDataDegX,eyeDataDegY] = convertEyeDataPix2DegML(eyeGazeData)
% given pixel values 
% and the following fixed variables:
% resolution of the monitor: 1280x720 pixel
% dimentions of the monitor: 532.5x300 cm
% change it to the degrees

xResolution = 1280;
yResolution = 720;
xDimMonitorBenQ = 532.55;
yDimMonitorBenQ = 300;
subjectDistance = 580;

%
xPixelsPerMM = xDimMonitorBenQ/xResolution;
yPixelsPerMM = yDimMonitorBenQ/yResolution;

% givenData
xPixelToConvert = eyeGazeData(:,1);
yPixelToConvert = eyeGazeData(:,2);

xPixelFromCenter = xPixelToConvert-(xResolution/2);
yPixelFromCenter = yPixelToConvert-(yResolution/2);

% calculations
eyeDataDegX = rad2deg(atan((xPixelFromCenter*xPixelsPerMM)/subjectDistance));
eyeDataDegY = -rad2deg(atan((yPixelFromCenter*yPixelsPerMM)/subjectDistance)); % minus to flip the Y axis
end