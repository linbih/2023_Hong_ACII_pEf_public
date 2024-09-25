
% pEf - EEG Analyses - MBBI - Step I | load BVA processed GA removed
% datasets

% Based on test_170604_EEG_34elec_Step_i and ii, mbbi_180420_EEG_Step_i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORTANT NOTE %%%

% Note that the following code is for subjects post May.31, i.e. Sub10
% onwards

%%% FIN %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG Pre-processing for subjects AFTER May.31 2018

% Step 1 | Import BVA processed GA removed dataset

% Step 2 | Load PD .mat file into EEGLAB dataset

% Step 3 | Fuse RT from .mat files into EEGLAB dataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------- %
% Code function  %
% -------------- %

% This code is EEG pre-processing for BVA exported, GA removed datasets

% -------------- %
% Code history   %
% -------------- %

% Created by LH, 180621
% Added run # exception for Sub19, 180624

% Updated by LH, 180715
% Added in comments regarding RT fusing from PD .mat to EEGLAB structure

% Updated by LH, 180822
% Commented out AREA -> DIAMETER conversion, as it's already been taken
% care of in updated PD_Step_i.m

% Aldo modified input option on pop_importeyetracker(), so that we're
% loading the pre-processed and z-scored PD traces

% Also changed dataset name from '_with_raw_pd' to '_with_clean_pd',
% '_with_pd_rt' will subsequently reflect the clean PD, not the raw PD
% apparently

% Updated by LH, 181117
% specified run number exception for Sub23, raw eeg file naming exception
% for Sub24

% Updated by LH, 181124
% Updated code w.r.t. updated function name of mbbi_file_mapping()

% Updated by LH, 190101
% specified run index exception for Sub35, raw eeg file naming exception
% for Sub35

% Updated by LH, 190107
% commented out run_num variable's exception for Sub15+16, as we haven't
% started analyzing those and this is to avoid confusion in future analysis

% Updated by LH, 190305
% added in run number exception for Sub15, as we'll start from Run2 instead
% of Run1


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
filepath = '/Users/linbihong/Dropbox/LIINC_2016_April_Onwards/A4_Simultaneous_Pupil_EEG_fMRI/Data_Analysis/A2_EEG_Analyses';

% If folder for individual subject doesn't exist, then create one.
if exist([filepath '/' sub_name '/Pre-processed dataset'],'dir') ~= 7
    mkdir(filepath,['/' sub_name '/Pre-processed dataset']);
    filepath = [filepath '/' sub_name '/Pre-processed dataset'];
else
    filepath = [filepath '/' sub_name '/Pre-processed dataset'];
end

% Run # exception for Sub15+16, added by LH, 180620
% Run # exception for Sub19, added by LH, 180624

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

% BVA processed EEG filename

bva_filename = cell(1,run_num);

for run_idx = 1:run_num
        
    % Below sub11 exception added by LH, 180610
    % Below sub24 exception added by LH, 181117
    
   % Run idx exception for Sub35, added by LH, 190101
   
    % Run # exception for Sub15, added by LH, 190305  
    
    if strcmp(sub_code,'sub11')
        bva_filename{run_idx} = [sub_name(3:6) '_' num2str(run_idx) '_GA_Remov_noF_21vols.vhdr'];
    
    elseif strcmp(sub_code,'sub15')
        bva_filename{run_idx} = [sub_name(1:6) '_' num2str(run_idx+1) '_GA_Remov_noF_21vols.vhdr'];
        
    elseif strcmp(sub_code,'sub19')
        bva_filename{run_idx} = [sub_name(1:6) '_' num2str(run_idx+1) '_GA_Remov_noF_21vols.vhdr'];
    
    elseif strcmp(sub_code,'sub24')
        bva_filename{run_idx} = [sub_name(1:6) '_2_' num2str(run_idx) '_GA_Remov_noF_21vols.vhdr'];
    
    elseif strcmp(sub_code,'sub35')
        bva_filename{run_idx} = [sub_name(1:6) '_2_' num2str(run_idx+1) '_GA_Remov_noF_21vols.vhdr'];
    
    else
        bva_filename{run_idx} = [sub_name(1:6) '_' num2str(run_idx) '_GA_Remov_noF_21vols.vhdr'];
    end
       
end

%%

% ------------------------------------------- %
%%% Load BVA processed GA-removed datasets %%%
% ------------------------------------------- %

for run_idx = 1:run_num
        
    [EEG,com] = pop_loadbv(filepath,bva_filename{run_idx});        
    
	eeglab_filename = [sub_code '_run' num2str(run_idx) '_bva_ga_removed'];
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Step 1 | Load %s to EEGLAB' ...
    '\n--------------------------------------------------------------\n'],eeglab_filename);
    
    EEG = pop_editset(EEG,'setname',eeglab_filename);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);
    eeglab r;

    % Save dataset
    EEG = pop_saveset(EEG,'filename',[eeglab_filename '.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
    %%% Following updated by LH, 180620
    
    % ------------------------------------------- %
    %%% Load PD .mat file %%%
    % ------------------------------------------- %
    
    % Import in PD dataset
    % Inputs - (EEG, .mat file to load, [start-event end-event],[which columns to
    % import], {new label names}, importEyeEvents, doRegression,
    % filterEyetrack, plotFig)
    
    % initial filepath
    filepath_init = '/Users/linbihong/Dropbox/LIINC_2016_April_Onwards/A4_Simultaneous_Pupil_EEG_fMRI/Data_Analysis';
    filepath_save = [filepath_init '/A1_PD_Analyses/' sub_name '/Pre-processed dataset/'];
    dataset_name_sti = [sub_code '_run' num2str(run_idx) '_pd_sti_clean_zs'];
    dataset_name_rt = [sub_code '_run' num2str(run_idx) '_pd_rt'];
    
    EEG = pop_importeyetracker(EEG,[filepath_save dataset_name_sti '.mat'], [8 8], 1:5, {'TIME' 'R-GAZE-X' 'R-GAZE-Y' 'R-DIAMETER' 'INPUT'}, 0, 1, 0, 0, 4);
    
    % Store the dataset into EEGLAB
    [ALLEEG EEG CURRENTSET ] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    
    % save into a new dataset
    EEG = pop_saveset(EEG,'filename',[eeglab_filename '_with_clean_pd.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
%     % convert area to diameter
%     % load channel 68, which corresponds to area
%     
%     area_label = EEG.chanlocs(68).labels;
%     if strcmp(area_label(end-3:end),'AREA')
%         EEG.data(68,:) = 256 * sqrt(EEG.data(68,:)/pi);
%     end
    
    % add RT event to the dataset
    
    pd_sti_mat = load([filepath_save dataset_name_sti '.mat']);
    pd_rt_mat = load([filepath_save dataset_name_rt '.mat']);
    
    % find first standard stimulus's latency [IN PD SAMPLES]
    first_std = pd_sti_mat.event(1,1);
    
    % find the 'latency' in samples (per EEGLAB's structure) to facilitate
    % RT event insertion in a bit
    
    % temp is first event of this type in ms
    [temp,all_sq_lat_std] = eeg_getepochevent(EEG,'S  8');
    
    % find first standard stimulus's latency [IN EEG SAMPLES]
    first_std_latency = EEG.srate * (temp/1000) + 1;
    
    % align RT events relative to the first standard stimulus 
    
    % actual time difference in secs - [IN PD
    % SAMPLES / PD SAMPLING RATE of 1KHz]
    rt_first_std = (pd_rt_mat.event(:,1) - first_std) / 1000;
    
    % insert RT into EEGLAB .set structure
    
    nevents = length(pd_rt_mat.event);
    
    for index = 1 : nevents
        
        % Add events relative to existing events
        
        % Add event to end of event list
        EEG.event(end+1) = EEG.event(index); 
        
        % Specifying the event latency [in EEG samples], to be whatever 
        % samples after the referent - the first standard 
        % stimulus - event. 
        
        EEG.event(end).latency = first_std_latency + rt_first_std(index)*EEG.srate; 
        EEG.event(end).type = 'Button'; 
    
    end

    % Check all events for consistency
    EEG = eeg_checkset(EEG, 'eventconsistency'); 
    
    % Store the dataset into EEGLAB
    [ALLEEG EEG CURRENTSET ] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    
    % save into a new dataset
    EEG = pop_saveset(EEG,'filename',[eeglab_filename '_with_pd_rt.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
    %%% Above updated by LH, 180620
    

end
