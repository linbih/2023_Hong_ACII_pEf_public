
% pEf - EEG Analyses - MBBI - Step IV.ica | Reject trials w/ trials identified after inspecting BOTH EEG + PD blink IC removed data; Concatenate over runs

% Based on mbbi_180621_EEG_Step_iv + 180927_EEG_Step_iv_v2.m

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

% --------------------------- %
% Dataset loading order       %
% --------------------------- % 

% Dataset 1 | sti-std dataset
%         2 | sti-odd dataset

% --------------------------------------------------------------- %
% Load GA removed + filtered + epoched datasets %
% --------------------------------------------------------------- %

sti_type = {'std','odd'};

for run_idx = 1:run_num

    for idx = 1:2

        close all;
        eeglab;
        
        dataset_name = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed sti-' sti_type{idx}];        
        eeg_filename = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed'];

        fprintf(['\n--------------------------------------------------------------' ...
        '\n   Load %s to EEGLAB' ...
        '\n--------------------------------------------------------------\n'],dataset_name);    

        EEG = pop_loadset([dataset_name '.set'],filepath);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, idx);
        eeglab redraw;

        % -------------------------- %
        % Reject trials %
        % -------------------------- %
        % 

        % ------------------------------------------ %
        % Apply semi-automatic epoch rejection      %
        % ------------------------------------------ %

        % Criterion: single channel - 6 s.d.; all channel - 2 s.d.; EEGPLOT type
        % Inputs: INEEG, typerej, elec_comp, locthresh, globthresh, superpose, reject, vistype

        % Manually inspect, then 'Update marks'; Since we also need to reject the incorrect response trials, 
        % use pause right after automatic rejection, so that when dataset is plotted again, we can see 
        % previously rejected trial info (before manually rejecting the incorrect response trials if needed)        
        
        EEG = pop_jointprob(EEG, 1, 1:63, 6, 2, 1, 0, 1, [], 1);
        pause();

        % Inputs: pop_eegplot( EEG, icacomp, superpose, reject )
        pop_eegplot(EEG, 1, 1, 0);
        pause();
        
        EEG = pop_editset(EEG,'setname',[eeg_filename ' sti-' sti_type{idx} '_epoch_marked']);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);      
       
        % we use reject.rejmanual instead of reject.rejjp here for we did
        % probability distribution based rejection FIRST, manual rejection
        % SECOND; hence manual rejection contains the final information
        
        % Also note that if there's no trial to reject, trial_rej_info
        % would be an empty matrix        
        
        trial_rej_info = EEG.reject.rejmanual;
        save(fullfile(filepath,[EEG.setname '.mat']),'trial_rej_info'); 
        
       % ------------------------------------- %
        % Load trial rejection info from EEG %
        % ------------------------------------- %
        %         
        
        % Extract indices of nonzero elements, i.e. trial indices to be
        % rejected
        
        trial_rej_mat = load([filepath '/' [eeg_filename ' sti-' sti_type{idx} '_epoch_marked.mat']]);
        trial_rej_info_eeg = find(trial_rej_mat.trial_rej_info);                

        % ------------------------------------- %
        % Load trial rejection info from PD %
        % ------------------------------------- %
        %                 
        
        trial_rej_info_pd_allrun = mbbi_trial_rej_pd_mapping(sub_code,sti_type{idx});

        trial_rej_info_pd = trial_rej_info_pd_allrun{run_idx};           
        
        trial_rej_all = union(trial_rej_info_eeg,trial_rej_info_pd);
        
        % if trial_rej_all is a column vector, we need to transpose it to a
        % row vector for the disp() to work (otherwise there'll be an
        % error) - commented by LH, 180928
        
        if size(trial_rej_all,2) == 1
            disp(['Rejected trial (PD+EEG): ',num2str(trial_rej_all')]);
        else
            disp(['Rejected trial (PD+EEG): ',num2str(trial_rej_all)]);
        end

        save(fullfile(filepath,[EEG.setname '_pd_eeg.mat']),'trial_rej_all'); 

        % -------------------------- %
        % Reject trials %
        % -------------------------- %
        %     
        
        % current dataset number should corresponds to idx + 1
        
        if ~strcmp(ALLEEG(idx).filename,[eeg_filename ' sti-' sti_type{idx} '.set'])
            current_dataset_num = str2double(input('Dataset number for epoched EEG dataset: ','s'));
        else
            current_dataset_num = idx;
        end

        EEG = pop_rejepoch(ALLEEG(current_dataset_num),trial_rej_all);

        EEG = pop_editset(EEG,'setname',[eeg_filename ' sti-' sti_type{idx} '_epoch_rejected_pd_eeg']);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);

        EEG = pop_saveset(EEG,'filename',[eeg_filename ' sti-' sti_type{idx} '_epoch_rejected_pd_eeg.set'],'filepath',filepath);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
        eeglab redraw;     
        
        % for oddball datasets - extract RT
        % Added by LH, 180623
        if idx == 2
            [temp,all_sq_lat_rt] = eeg_getepochevent(EEG,'Button');
            
            % RT in seconds
            rt = cell2mat(all_sq_lat_rt) / 1000;
            save(fullfile(filepath,[EEG.setname '_RT.mat']),'rt'); 
        end                         

    end

end

%% Concatenate all runs

for idx = 1:2

    close all;
    eeglab;
    
    run_num_start = 1;
    
    for run_idx = run_num_start:run_num

        dataset_name = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed sti-' sti_type{idx} '_epoch_rejected_pd_eeg'];
        
        fprintf(['\n--------------------------------------------------------------' ...
        '\n   Load %s to EEGLAB' ...
        '\n--------------------------------------------------------------\n'],dataset_name);    

        EEG = pop_loadset([dataset_name '.set'],filepath);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, run_idx);

    end

    eeglab redraw;
    
    % use run_num_start to include all eligible runs together, 180717
    EEG= pop_mergeset(ALLEEG,run_num_start:run_num);
    
    % added exception clause for Sub15 - LH, 180717
    % added exception clause for Sub23 - LH, 181120
    
    % optimized Run # exception clause expression (without specifying
    % subjects - as that's already been taken care of by run_num variable),
    % LH, 190107    
    
    all_run_filename = [sub_code '_' num2str(run_num)  ...
                        'runs_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed sti-' ...
                        sti_type{idx} '_epoch_rejected_pd_eeg'];        
    
    EEG = pop_editset(EEG,'setname',all_run_filename);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);
    eeglab redraw;

    EEG = pop_saveset(EEG,'filename',[all_run_filename '.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
end

