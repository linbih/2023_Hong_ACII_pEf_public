### This repository contains code for:

Hong L, He H, Sajda P (2023). Capturing Interactions between Arousal and Cortical Dynamics with Simultaneous Pupillometry and EEG-fMRI. 2023 11th International Conference on Affective Computing and Intelligent Interaction (ACII). [Link](https://ieeexplore.ieee.org/abstract/document/10388080)


### Data can be accessed at:
- gradient artifact removed EEG (and pupil) data at: [Google Drive Loc1](https://drive.google.com/file/d/1vPaoCBM-3AyWuON679MZzqAMCDzDiONG/view?usp=sharing)
- pre-processed (BCG artifact removed, re-referenced, ICA-based artifacts removed) EEG data at: [Google Drive Loc2](https://drive.google.com/file/d/1mA_rSV7Q_fl1nrfxsEjT7hZq2FEiB0mh/view?usp=sharing)
- fMRI data at: [Google Drive Loc3](https://drive.google.com/drive/folders/1rK1oeFc7imxYXffMCfQ3QBPQ_qZQkk4f?usp=sharing)

### Example use case: 
1. perform single trial EEG analysis
```
% as seen on line 146 of mbbi_EEG_Step_v_ica.m

% run logistic regression and save results in mat file
[ALLEEG, LR] = runLR(ALLEEG, setlist, chansubset, ...
    subjectID, condition, ...
    domainname, windowlength, windowoffset, LOO, filepath_lr, epochtype);

% plot the results Az vs time
plotAz(LR,[],[],Azplotdir,1);
```
  
### More methodological details at:
Hong L, He H, Sajda P (2022). Pupil-linked phasic arousal relates to the reduction of response inhibition: inferences from a simultaneous pupillometry and EEG-fMRI study. Pre-print available on [Biorxiv](https://www.biorxiv.org/content/10.1101/2022.08.22.504728v1).

