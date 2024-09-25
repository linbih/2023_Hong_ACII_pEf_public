
% pEf - EEG Analyses - MBBI - Step VII.ica | Generate single-trial EV files w/ trials identified after inspecting BOTH EEG + PD ICA-processed data;

% Based on mbbi_190325_EEG_Step_vi_ica.m

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

% Created by LH, 190409

% Updated by LH, 190410
% options on regressor height

% Updated by LH, 190522
% ran in a loop for all subs, using ori option, though commented out afterwards

%%

% ------------------------------------------- %
%%% Specify environment %%%
% ------------------------------------------- %

close all;
eeglab
      
%[sub_code_all, group_code] = mbbi_group_code_setting();

%for idx = 1:length(sub_code_all)

%sub_code = sub_code_all{idx};
sub_code = input('Which subject will we be looking at? ','s');

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


% ------------------------------------------- %
%%% Load RT + Y-value %%%
% ------------------------------------------- %

%%% RT

rt = [];

for run_idx = 1:run_num

    dataset_name = [sub_code '_run' num2str(run_idx) '_withQRS_BCGremoved_sansECGchan_re-ref_withICA_blinkIC_removed'];

    rt_mat = load([filepath '/' [dataset_name ' sti-odd_epoch_rejected_pd_eeg_RT.mat']]);
    rt_info = rt_mat.rt;

    rt = [rt;rt_info'];

end

%%% y-value

eeg_filename = ['LRresults_' sub_name '__Sti-locked_Odd_vs_Std_|_BCG_artifact_+_blink_IC_removed_|_PD+EEG_Auditory.mat'];

y_value = load([filepath '/LR_results/' eeg_filename]);

%%% why LR.yLOO{1} is always odd

% input for runLR()                     | setlist(1) - odd; setlist(2) - std
% input for jen_logisticregression()    | setlist(1) - std; setlist(2) - odd / due to use of fliplr()
% output for jen_logisticregression()   | yLOO(1) - std; yLOO(2) - odd
% output for runLR()                    | LR.yLOO(1) - odd; LR.yLOO(2) - std / due to line 256-59

y_odd = y_value.LR.yLOO{1};
y_std = y_value.LR.yLOO{2};

%%% what y-value to plot
% 1 - y-value as is
%     | {y_std,y_odd}

% 2 - de-meaned y-value (mean of all trials per class per certain window)
%     | {y_std_demean,y_odd_demean}

% 3 - y-value w/ RT regressed out then de-meaned
%     | {y_std_demean,res_y_odd_demean}

% 4 - z-scored y-value (z-scored per class per certain window)
%     | {y_std_zs,y_odd_zs}

% 5 - y-value w/ RT regressed out then z-scored 
%     | {y_std_zs,res_y_odd_zs}

y_value_idx = input('What type of y-value to modulate regressor height (1/2/3/4/5 | ori/ori_demean/res_demean/ori_zs/res_zs) ','s');
y_value_filesave_type = {'ori','ori_demean','res_demean','ori_zs','res_zs'};

if strcmp(y_value_idx,'1')
    
    % 1 - y-value as is
   
    amplitude_y = {y_std,y_odd};
    
elseif strcmp(y_value_idx,'2')
    
    % 2 - de-meaned y-value (mean of all trials per class per certain window)
    
    y_odd_demean = zeros(size(y_odd,1),size(y_odd,2));
    y_std_demean = zeros(size(y_std,1),size(y_std,2));
    
    for winidx = 1:size(y_odd,2)  
        mean_y_odd = mean(y_odd);
        mean_y_std = mean(y_std);
        
        y_odd_demean(:,winidx) = y_odd(:,winidx) - mean_y_odd(winidx);
        y_std_demean(:,winidx) = y_std(:,winidx) - mean_y_std(winidx);
    end
    
    amplitude_y = {y_std_demean,y_odd_demean};
    
elseif strcmp(y_value_idx,'3')
    
    % 3 - y-value w/ RT regressed out then de-meaned    
    
    res_y_odd = zeros(size(y_odd,1),size(y_odd,2));

    for window_idx = 1:size(y_odd,2)
        p = polyfit(rt,y_odd(:,window_idx),1);
        yfit = polyval(p,rt);

        res_y_odd(:,window_idx) = y_odd(:,window_idx) - yfit;
    end  

    res_y_odd_demean = zeros(size(res_y_odd,1),size(res_y_odd,2));
    y_std_demean = zeros(size(y_std,1),size(y_std,2));
    
    for winidx = 1:size(y_odd,2)  
        mean_y_odd = mean(res_y_odd);
        mean_y_std = mean(y_std);
        
        res_y_odd_demean(:,winidx) = res_y_odd(:,winidx) - mean_y_odd(winidx);
        y_std_demean(:,winidx) = y_std(:,winidx) - mean_y_std(winidx);
    end
    
    amplitude_y = {y_std_demean,res_y_odd_demean};
    
elseif strcmp(y_value_idx,'4')
    
    % 4 - z-scored y-value (z-scored per class per certain window)
    
    y_odd_zs = zscore(y_odd);
    y_std_zs = zscore(y_std);
    
    amplitude_y = {y_std_zs,y_odd_zs};   
    
elseif strcmp(y_value_idx,'5')
    
    % 5 - y-value w/ RT regressed out then z-scored

    y_std_zs = zscore(y_std);
    
    res_y_odd = zeros(size(y_odd,1),size(y_odd,2));

    for window_idx = 1:size(y_odd,2)
        p = polyfit(rt,y_odd(:,window_idx),1);
        yfit = polyval(p,rt);

        res_y_odd(:,window_idx) = y_odd(:,window_idx) - yfit;
    end  
    
    res_y_odd_zs = zscore(res_y_odd);
    
    amplitude_y = {y_std_zs,res_y_odd_zs};    
    
end  


%%

sti_type = {'std','odd'};

std_trial_num = zeros(run_num_start,run_num+1);
odd_trial_num = zeros(run_num_start,run_num+1);

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
    %%% STV fMRI Analyses EV files %%%
    % ------------------------------------------------------------ %     

    % EV files for both standard and oddball

    for idx = 1:2

        % ------------------------------------------- %
        %%% Tease out trials already rejected      %%%
        % ------------------------------------------- %    

        trial_rej_mat = load([filepath '/' [dataset_name ' sti-' sti_type{idx} '_epoch_marked_pd_eeg.mat']]);
        trial_rej_info = trial_rej_mat.trial_rej_all;
        
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
            
            std_trial_num(run_idx+1) = std_trial_num(run_idx) + length(std_latency_new);

        else
            
            if ~isempty(trial_rej_info)
                odd_latency_new = odd_latency;
                odd_latency_new(trial_rej_info) = [];
            else
                odd_latency_new = odd_latency;
            end
            
            odd_trial_num(run_idx+1) = length(odd_latency_new);        

            odd_trial_num(run_idx+1) = odd_trial_num(run_idx) + length(odd_latency_new);            
            
        end
                
        % ------------------------------------------- %
        %%% Generate files      %%%
        % ------------------------------------------- %  
        
        for winidx = 1:length(y_value.LR.windowcenter)
                
            if idx == 1
                y_value_type = {'ori','ori_demean','ori_demean','ori_zs','ori_zs'};
            else
                y_value_type = {'ori','ori_demean','res_demean','ori_zs','res_zs'};
            end
            
            fname = [sub_code '_run' num2str(run_idx) ...
                     '_BCGremoved_blinkIC_removed_rejected_pd_eeg_EV_STV_' ...
                     sti_type{idx} '_y_' y_value_type{str2num(y_value_idx)} '_' num2str(y_value.LR.windowcenter(winidx)) '_ms'];
                        
            % save path to our directory of choice

            imgfilepath = [filepath_init '/' sub_name '/STV_files/' y_value_filesave_type{str2num(y_value_idx)}];

            % If folder doesn't exist, then create one.
            if exist(imgfilepath,'dir') ~= 7
                mkdir(imgfilepath);
            end

            fullfname=[imgfilepath '/' fname '.txt'];
            fid = fopen(fullfname,'w');        

            %%% Stimulus onset - 
            %%%                     std-latency + windowcenter - 1/2 * 100ms for sti-std, 
            %%%                     odd-latency + windowcenter - 1/2 * 100ms for sti-odd, 

            %%% Duration - 0.1 for STV EVs
            duration = 0.1;
            
            if idx == 1
                stimulus_latency = std_latency_new + y_value.LR.windowcenter(winidx)/1000 - 0.5 * duration;
            else
                stimulus_latency = odd_latency_new + y_value.LR.windowcenter(winidx)/1000 - 0.5 * duration;
            end
            
            %%% Amplitude - y-value

            if idx == 1
                trial_num = std_trial_num;
            else
                trial_num = odd_trial_num;
            end
            
            trial_idx = (trial_num(run_idx)+1):trial_num(run_idx+1); 
            amplitude = amplitude_y{idx}(trial_idx,winidx);
            
            fprintf(['\n--------------------------------------------------------------' ...
            '\n   Stimulus length: %d' ...
            '\n--------------------------------------------------------------\n'],length(stimulus_latency));        
        
            % If indeed # of trials match between y-value and stimulus
            % variable (double check)
            
            if length(amplitude) == length(stimulus_latency)                      
                
                for i = 1:length(stimulus_latency)
                    % write to file
                    fprintf(fid(1),'%6.3f %6.3f %6.3f',[stimulus_latency(i) duration amplitude(i)]);

                    if i ~= length(stimulus_latency)
                        fprintf(fid,'\n');
                    end
                end
                
            else
                fprintf(['\n--------------------------------------------------------------' ...
                '\n   Stimulus length DOES NOT MATCHES y-value length !! ' ...
                '\n--------------------------------------------------------------\n']);  

            end
                
            fclose(fid);

            clear fullfname
            clear fid
        
        end

    end
    
    close all    
    
end

% end
