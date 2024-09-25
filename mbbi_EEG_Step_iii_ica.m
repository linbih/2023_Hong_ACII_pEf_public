
% pEf - EEG Analyses - MBBI - Step III.ica | Remove artifactual IC from all-run datasets (w/ ICA results) + divide all-run dataset into single-run dataset + epoch

% Based on mbbi_180621_EEG_Step_iii.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORTANT NOTE %%%

% Note that the following code is for subjects post May.31, i.e. Sub10
% onwards

%%% FIN %%%

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

% Updated by LH, 190322
% added 'res-flux', i.e. ERP check for each run

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

if exist([filepath_init '/ICA_Analyses'],'dir') ~= 7
    disp('Directory to load dataset doesn''t exist!');
%     break;
else
    filepath_load = [filepath_init '/ICA_Analyses'];
end

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

% ------------------------------------------- %
%%% Load datasets w/ ICA results %%%
% ------------------------------------------- % 

dataset_name = [sub_code '_' num2str(run_num) 'runs_withQRS_BCGremoved_sansECGchan_re-ref_withICA'];

EEG = pop_loadset('filename',[dataset_name '.set'], ...
                  'filepath',filepath_load);

EEG = eeg_checkset(EEG);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;

% ------------------------------------------- %
%%% Reject IC of eye blink %%%
% ------------------------------------------- % 

% extract blink IC for each subject
blink_ic_idx = blink_ic_mapping(sub_code);

% remove blink IC
EEG = pop_subcomp(EEG, blink_ic_idx, 0);
EEG = pop_editset(EEG,'setname',[dataset_name '_blinkIC_removed']);
[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);
eeglab redraw;

current_dataset_num = 2;

if strcmp(ALLEEG(current_dataset_num).setname,[dataset_name '_blinkIC_removed'])
    eeg_idx = current_dataset_num;
else
    eeg_idx = str2double(input('Dataset number for blink IC removed dataset: ','s'));
end  

%%

% ------------------------------------------- %
%%% Chop all-run dataset into single-run %%%
% ------------------------------------------- % 

allEventTypes = {ALLEEG(eeg_idx).event.type}';
allEventCodes = {ALLEEG(eeg_idx).event.code}';

keywordIdx = intersect(find(strcmp(allEventTypes, 'boundary')),find(strcmp(allEventCodes, 'New Segment')));

if length(keywordIdx) == run_num        
    
    for run_idx =  1:run_num
        
        if run_idx ~= run_num
            run_dur = [ALLEEG(eeg_idx).event(keywordIdx(run_idx)).latency ...
                       ALLEEG(eeg_idx).event(keywordIdx(run_idx+1)).latency - 1];
                   
        else
            run_dur = [ALLEEG(eeg_idx).event(keywordIdx(run_idx)).latency ...
                       ALLEEG(eeg_idx).pnts];    
        end
    
        EEG = pop_select(ALLEEG(eeg_idx),'point',run_dur);
        
        EEG = pop_editset(EEG,'setname',[sub_code '_run' num2str(run_idx) ...
                                         '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed']);
        
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);
        eeglab redraw; 

        EEG = pop_saveset(EEG,'filename',[sub_code '_run' num2str(run_idx) ...
                                         '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed'], ...
                              'filepath',filepath);
        
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
        eeglab redraw;        
        
        % Below added by LH, 190322
        % compare ERP results w/ previous approach
        
        filepath_plot = '/Users/linbihong/Dropbox/LIINC_2016_April_Onwards/A4_Simultaneous_Pupil_EEG_fMRI/Data_Analysis/A2_EEG_Analyses';

        % If folder for individual subject doesn't exist, take notice
        if exist([filepath_plot '/' sub_name '/Pre-processed dataset'],'dir') ~= 7
            disp('Directory to load dataset doesn''t exist!');
            %break;
        else
            filepath_plot = [filepath_plot '/' sub_name '/Pre-processed dataset'];
        end
        
        mbbi_res_flux_est_ica(ALLEEG, EEG, (eeg_idx+run_idx), filepath_plot, [sub_code '_run' num2str(run_idx)]);

        
    end
    
end

%%

% ------------------------------------------- %
%%% Carry out epoching %%%
% ------------------------------------------- % 

for run_idx = 1:run_num

    close all;
    eeglab;
    
    % ------------------------------------------- %
    %%% Load single-run blink IC removed dataset %%%
    % ------------------------------------------- %    
    
    dataset_name = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed'];
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Load %s to EEGLAB' ...
    '\n--------------------------------------------------------------\n'],dataset_name);    
    
    EEG = pop_loadset([dataset_name '.set'],filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;    
    
    % ------------------ %
    % Epoching  %
    % ------------------ %        
    
    current_dataset_num = 1;

    if strcmp(ALLEEG(current_dataset_num).filename,[dataset_name '.set'])
        eeg_idx = current_dataset_num;
    else
        eeg_idx = str2double(input('Dataset number for blink IC removed EEG dataset: ','s'));
    end      
    
    % Epoch Stimulus-locked standard trials
    EEG = pop_epoch(ALLEEG(eeg_idx),{'S  8'},[-0.5 2],'newname',[dataset_name ' sti-std']);
    EEG = pop_rmbase(EEG,[-500 0]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);

    EEG = pop_saveset(EEG,'filename',[EEG.setname '.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
    % Epoch Stimulus-locked oddball trials
    EEG = pop_epoch(ALLEEG(eeg_idx),{'S 32'},[-0.5 2],'newname',[dataset_name ' sti-odd']);
    EEG = pop_rmbase(EEG,[-500 0]);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);

    EEG = pop_saveset(EEG,'filename',[EEG.setname '.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;    

end
