
% pEf - EEG Analyses - MBBI - Step V.ica | Run LR!! on blink IC removed PD+EEG combined criteria-based trial rejected dataset

% Based on mbbi_180927_EEG_Step_v_v2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORTANT NOTE %%%

% Note that the following code is for subjects post May.31, i.e. Sub10
% onwards

%%% FIN %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG Processing for subjects AFTER May.31 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------- %
% Code function  %
% -------------- %

% This code is EEG Processing for BVA exported, GA removed + properly
% filtered datasets

% -------------- %
% Code history   %
% -------------- %

% Created by LH, 190321

%%

% ------------------------------------------- %
%%% Specify environment %%%
% ------------------------------------------- %

close all;
eeglab
                     
sub_code = input('Which subject will we be looking at? ','s');

% Go into corresponding folder for specified subject
sub_name = mbbi_file_mapping(sub_code);

% initial filepath
filepath_init = '/Users/linbihong/Dropbox/LIINC_2016_April_Onwards/A4_Simultaneous_Pupil_EEG_fMRI/Data_Analysis/A2_EEG_Analyses';

% If folder for individual subject doesn't exist, then create one.
if exist([filepath_init '/' sub_name '/STV_dataset'],'dir') ~= 7
    disp('Directory to load dataset doesn''t exist!');

else
    filepath = [filepath_init '/' sub_name '/STV_dataset'];
end

% Run # exception for Sub15+16
% Run # exception for Sub23, added by LH, 181117

% Run # exception for Sub15+16 commented out, LH, 190107

% Run # exception for Sub15, added by LH, 190305

% if strcmp(sub_code,'sub15') || strcmp(sub_code,'sub16') || strcmp(sub_code,'sub19')
if strcmp(sub_code,'sub19')
    run_num = 4;
elseif strcmp(sub_code,'sub15') || strcmp(sub_code,'sub23')    
    run_num = 3;
    
else
    run_num = 5;
end

%%

% --------------------------------------------------------------- %
% Load GA removed + filtered + epoched + trial rejected datasets %
% --------------------------------------------------------------- %

sti_type = {'std','odd'};

for idx = 1:2

    % added exception clause for Sub23 - LH, 181120
	
    % optimized Run # exception clause expression (without specifying
    % subjects - as that's already been taken care of by run_num variable),
    % LH, 190107 
    
    dataset_name = [sub_code '_' num2str(run_num) ...
                    'runs_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed sti-' ...
                    sti_type{idx} '_epoch_rejected_pd_eeg'];        
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Step 8.0 | Load %s to EEGLAB' ...
    '\n--------------------------------------------------------------\n'],dataset_name);    

    EEG = pop_loadset([dataset_name '.set'],filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, idx);

end

eeglab redraw;

%%

% --------------------------------------------------- %
% Parameter configuration for LR  %
% --------------------------------------------------- %

domainname = 'time';    % or 'frequency' (for logistic regression)

% setlist = [setlist_odd setlist_std]
% and here we have setlist_odd as dataset 2; setlist_std as dataset 1
setlist = [2 1];  
chansubset = 1:length(EEG.chanlocs); % channels to use for LR

subjectID = sub_name;
condition = 'Sti-locked Odd vs Std | BCG artifact + blink IC removed | PD+EEG';

windowlength = 50; % width of training window
windowoffset = 0:25:1000; % centers of training window relative to event

LOO = 1;                % perform leave-one-out cross validation
filepath_lr = filepath;

% updated variable name from sessiontype to epochtype, LH, 190108
epochtype = 'Auditory';
Azplotdir = [filepath '/Az_plot'];

fprintf(['\n--------------------------------------------------------------------------------------' ...
    '\n   Training Classifier on' ...
    '\n   %s Data' ...
    '\n   (in %s domain)' ...
    '\n   Set names are %s ' ...
    '\n             and %s' ...
    '\n--------------------------------------------------------------------------------------\n'], ...
    dataset_name, domainname, ALLEEG(1).setname, ALLEEG(2).setname)

% --------------------------------------------------- %
% Run LR  %
% --------------------------------------------------- %

% run logistic regression and save results in mat file

[ALLEEG, LR] = runLR(ALLEEG, setlist, chansubset, ...
    subjectID, condition, ...
    domainname, windowlength, windowoffset, LOO, filepath_lr, epochtype);

% --------------------------------------------------- %
% Plot LR Results  %
% --------------------------------------------------- %

% plot the results Az vs time
plotAz(LR,[],[],Azplotdir,1);

% close figures and clear the EEGLAB datasets for next session type
close all
drawnow
ALLEEG = pop_delset(ALLEEG, 1:length(ALLEEG));
