This folder contains programs that were used to segment the data and to find bad trials and electrodes.

Codes written below are provided only for completeness. The eventual outcome of the codes is the generation of data that is available in the segmentedData folder. Since this data is already provided, these codes do not need to be run again

1. First, rawData must be saved in the data folder. The code saveRawData was used to do it.
2. Start by running runSegmentAndSaveData, which uses sagmentAndSaveData for doing the segmentation. Note that this code uses programs in commonPrograms available at https://github.com/supratimray/CommonPrograms, EEGLab as well as other dependencies to deal with BrainProducts data file.
3. Data segmentation failed for some specific reasons for a few cases (usually because markers were unavailable). Details are provided in the runSegmentAndSaveData file itself.
4. Next, run runFindBadTrialsWithEEG to use our automated preprocessing code to find a set of bad electrodes and a set of bad stimulus repeats that are common for all electrodes. Programs for this pipeline can be found in commonPrograms. The logic is explained in Murty and Ray, 2022, bioProtocol. The output of this program is a mat file called badTrials_wo_v8 within the segmentedData folder within each subject folder.

In addition, the program displayMeditationData in commonAnalysisCodes/displayCodes can be used to display preliminary results for a single subject. The program runDisplayMeditationData can be used to save all these results as .tif files, which are saved locally in this folder.

