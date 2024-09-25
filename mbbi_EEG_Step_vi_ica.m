
% pEf - EEG Analyses - MBBI - Step VI.ica | Generate EV files w/ trials identified after inspecting BOTH EEG + PD ICA-processed data;

% Based on mbbi_181107_EEG_Step_vi_v2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORTANT NOTE %%%

% Note that the following code is for subjects post May.31, i.e. Sub10
% onwards

%%% FIN %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------- %
% Code function  %
% -------------- %

% This code is EV file generation

% -------------- %
% Code history   %
% -------------- %

% Created by LH, 190325

% Updated by LH, 20'0108
% created another option of stimulus EV duration of 100ms, also run this in
% a loop for all subs

%%

% call function for group code setting

[sub_code_all, group_code] = mbbi_group_code_setting();

for sub_idx = 1:length(sub_code_all)

    % ------------------------------------------- %
    %%% Specify environment %%%
    % ------------------------------------------- %

    close all;
    eeglab

    % following line added by LH, 20'0108
    ev_duration = '0.1'; %input('Set the stimulus EV duration to 0.1 or 0.2 second? ', 's');

    sub_code = sub_code_all{sub_idx}; %input('Which subject will we be looking at? ','s');

    % Go into corresponding folder for specified subject
    sub_name = mbbi_file_mapping(sub_code);

    % initial filepath
    filepath_init = '/Users/linbihong/Dropbox/LIINC_2016_April_Onwards/A4_Simultaneous_Pupil_EEG_fMRI/Data_Analysis/A2_EEG_Analyses';

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

    run_num_start = 1;


    %%

    sti_type = {'std','odd','rt'};

    for run_idx = run_num_start:run_num

        % ------------------------------------------- %
        %%% Load GA removed filtered re-ref datasets %%%
        % ------------------------------------------- %    

        dataset_name = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed'];

        fprintf(['\n--------------------------------------------------------------' ...
        '\n   Load %s to EEGLAB' ...
        '\n--------------------------------------------------------------\n'],dataset_name);    

        EEG = pop_loadset([dataset_name '.set'],filepath);
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        eeglab redraw;

        % ------------------------------------------- %
        %%% Find 1st TR, i.e. start of experiment %%%
        % ------------------------------------------- %    

        [start_of_exp,all_sq_lat_tr] = eeg_getepochevent(EEG,'S128');

        % ------------------------------------------- %
        %%% Find stimuli's timing                   %%%
        % ------------------------------------------- %    

        [temp,all_sq_lat_std] = eeg_getepochevent(EEG,'S  8');
        [temp,all_sq_lat_odd] = eeg_getepochevent(EEG,'S 32');

        std_latency = (cell2mat(all_sq_lat_std) - start_of_exp) / 1000;
        odd_latency = (cell2mat(all_sq_lat_odd) - start_of_exp) / 1000;


        % ------------------------------------------------------------ %
        %%% TRADITIONAL fMRI Analyses EV files %%%
        % ------------------------------------------------------------ %     

        % EV files for both standard and oddball

        % Have other options of setting the stimulus EV duration to 0.1s, but
        % don't touch RT EV, as there's no need to change the EV duration of
        % that

        % Updated by LH, 20'0108         

        for idx = 1:2
    %     for idx = 1:3

            % ------------------------------------------- %
            %%% Tease out trials already rejected      %%%
            % ------------------------------------------- %    

            if idx ~=3

                trial_rej_mat = load([filepath '/' [dataset_name ' sti-' sti_type{idx} '_epoch_marked_pd_eeg.mat']]);
                trial_rej_info = trial_rej_mat.trial_rej_all;

            % For RT
            else
                trial_rej_mat = load([filepath '/' [dataset_name ' sti-odd_epoch_marked_pd_eeg.mat']]);
                trial_rej_info = trial_rej_mat.trial_rej_all;

                rt_mat = load([filepath '/' [dataset_name ' sti-odd_epoch_rejected_pd_eeg_RT.mat']]);
                rt_info = rt_mat.rt;
            end

            % exclude indices of trials to include: since trial_rej_info marks
            % the indices of trials to be excluded, we set the corresponding element to []

            % If there's no trial rejected, then simply take all stimuli's
            % timing and feed into the EV file

            if idx == 1

                if ~isempty(trial_rej_info)
                    std_latency_new = std_latency;
                    std_latency_new(trial_rej_info) = [];
                else
                    std_latency_new = std_latency;
                end
            else

                if ~isempty(trial_rej_info)
                    odd_latency_new = odd_latency;
                    odd_latency_new(trial_rej_info) = [];
                else
                    odd_latency_new = odd_latency;
                end
            end


            % ------------------------------------------- %
            %%% Generate files      %%%
            % ------------------------------------------- %  

            fname = [sub_code '_run' num2str(run_idx) '_BCGremoved_blinkIC_removed_rejected_pd_eeg_EV_traditional_' sti_type{idx}];

            % save path to our directory of choice

            imgfilepath = [filepath_init '/' sub_name '/STV_files'];

            % If folder doesn't exist, then create one.
            if exist(imgfilepath,'dir') ~= 7
                mkdir(imgfilepath);
            end

            if idx ~= 3

                % Have other options of setting the stimulus EV duration to 0.1s
                % Updated by LH, 20'0108        

                fullfname=[imgfilepath '/' fname '_' ev_duration 's.txt'];
                fid = fopen(fullfname,'w');

            else

                fullfname=[imgfilepath '/' fname '.txt'];
                fid = fopen(fullfname,'w');

            end

            %%% Stimulus onset - std-latency for sti-std, odd-latency for
            %%% sti-odd

            if idx == 1
                stimulus_latency = std_latency_new;
            else
                stimulus_latency = odd_latency_new;
            end

            %%% Duration - 0.2 for stimulus, or RT for rt

            % Have other options of setting the stimulus EV duration to 0.1s
            % Updated by LH, 20'0108

            if idx ~=3
                duration = repmat(str2double(ev_duration),length(stimulus_latency),1);
                % duration = repmat(0.2,length(stimulus_latency),1);
            else
                duration = rt_info;
            end

            %%% Amplitude - unit amp of 1

            amplitude = 1;

            fprintf(['\n--------------------------------------------------------------' ...
            '\n   Stimulus length: %d' ...
            '\n--------------------------------------------------------------\n'],length(stimulus_latency));        

            for i = 1:length(stimulus_latency)
                % write to file
                fprintf(fid(1),'%6.3f %6.3f %6.3f',[stimulus_latency(i) duration(i) amplitude]);

                if i ~= length(stimulus_latency)
                    fprintf(fid,'\n');
                end
            end

            fclose(fid);

            clear fullfname
            clear fid

        end

        close all    

    end

end

