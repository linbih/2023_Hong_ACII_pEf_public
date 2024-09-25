
% pEf - EEG Analyses - MBBI - Step II | pre-processing RT inserted BVA GA removed datasets | BPF @0.5-50Hz

% Based on test_170604_EEG_34elec_Step_ii, mbbi_180421_EEG_Step_ii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% IMPORTANT NOTE %%%

% Note that the following code is for subjects post May.31, i.e. Sub10
% onwards

%%% FIN %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EEG Pre-processing for subjects AFTER May.31 2018

% Step 4 | Load RT inserted BVA processed GA removed datasets
%          Output: eeglab struct

% Step 5 | Extract first 64 channels and create dataset w/ EEG data only
%          (last five channels are PD)

% Step 6 | Apply median filter

% Step 7 | Apply bandpass filter
%          Output: [sub_code '_run' num2str(run_idx)
%          '_bva_ga_removed_with_rt_medF_bandpassF'].set

% Step 8 | Check res-flux 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------- %
% Code function  %
% -------------- %

% This code is EEG pre-processing for RT inserted, BVA exported, GA removed datasets

% -------------- %
% Code history   %
% -------------- %

% Created by LH, 180621
% Added run # exception for Sub19, 180624

% Updated by LH, 181117

% specified run number exception for Sub23

% Updated by LH, 181124
% Updated code w.r.t. updated function name of mbbi_file_mapping()

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

% If folder for individual subject doesn't exist, take notice
if exist([filepath '/' sub_name '/Pre-processed dataset'],'dir') ~= 7
    disp('Directory to load dataset doesn''t exist!');
    %break;
else
    filepath = [filepath '/' sub_name '/Pre-processed dataset'];
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

% Specify specs for Bandpass Filter

fcutlow = 0.5;   %low cut frequency in Hz
order = 4;    
fcuthigh = 50;   %high cut frequency in Hz

%%

% ------------------------------------------- %
%%% Load BVA processed GA-removed datasets %%%
% ------------------------------------------- %

for run_idx = 1:run_num
    
    close all;
    eeglab;
    
    % ------------------------------------------- %
    %%% Load EEGLAB datasets %%%
    % ------------------------------------------- %
    
    dataset_name = [sub_code '_run' num2str(run_idx) '_bva_ga_removed_with_pd_rt.set'];
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Load %s to EEGLAB' ...
    '\n--------------------------------------------------------------\n'],dataset_name);

    % if such file exists (if exists, output is 2; otherwise will be 0)
    
    if exist(fullfile(filepath,dataset_name), 'file') == 2
        EEG = pop_loadset(dataset_name,filepath);
    else
        disp('EEGLAB dataset doesn''t exist!');
        break;
    end
    
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;  

    % ------------------------------------------------ %
    %%% Extract EEG data from 69 (EEG+PD) channels %%%
    % ------------------------------------------------ %       
    
    EEG = pop_select(ALLEEG(1), 'nochannel', 65:69);
    EEG = eeg_checkset(EEG);

    EEG = pop_editset(EEG,'setname',[dataset_name(1:end-15) '_with_rt']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);      

    EEG = pop_saveset(EEG,'filename',[dataset_name(1:end-15) '_with_rt.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;    
    
    % ------------------------------------------------ %
    %%% Apply median filter to get rid of residue GA %%%
    % ------------------------------------------------ %   
    
    if strcmp(ALLEEG(2).filename,[dataset_name(1:end-15) '_with_rt.set'])
        dataset = 2;
    else
        dataset = str2double(input('with-RT EEG dataset to look at: ','s'));
    end    
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Applying median filter on %s' ...
    '\n--------------------------------------------------------------\n'],dataset_name);   
    
    mbbi_data_medf = zeros(size(ALLEEG(dataset).data,1),size(ALLEEG(dataset).data,2));

    for chan = 1:size(ALLEEG(dataset).data,1)
        mbbi_data = ALLEEG(dataset).data(chan,:);

        mbbi_data_medf(chan,:) = medfilt1(double(mbbi_data),10);
        fprintf('Filtering chan %d\n',chan);
    end

    EEG.data = mbbi_data_medf;  
    
    % ------------------------------------------------ %
    %%% Apply bandpass filter to further clean the data %%%
    % ------------------------------------------------ %      
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Apply bandpass filter ' ...
    '\n--------------------------------------------------------------\n']);

    [b,a] = butter(order,[fcutlow,fcuthigh]/(EEG.srate/2),'bandpass');
    
    % figure(); 
    % freqz(b,a,1000,EEG.srate);

    EEG.data = filtfilt(b,a,double(EEG.data'))';

    % ---------------------------- %
    %%% Save filtered data %%%
    % ---------------------------- %         
            
    EEG = pop_editset(EEG,'setname',[dataset_name(1:end-15) '_with_rt_medF_bandpassF']);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG);      
    
    EEG = pop_saveset(EEG,'filename',[dataset_name(1:end-15) '_with_rt_medF_bandpassF.set'],'filepath',filepath);
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG,EEG,CURRENTSET);
    eeglab redraw;
    
    fprintf(['\n--------------------------------------------------------------' ...
    '\n   Bandpass filtered dataset saved in EEGLAB: %s ' ...
    '\n--------------------------------------------------------------\n'],EEG.setname);      

    % ---------------------------- %
    %%% Check res-flux for each run %%%
    % ---------------------------- %         
            
    if strcmp(ALLEEG(3).filename,[dataset_name(1:end-15) '_with_rt_medF_bandpassF.set'])
        dataset = 3;
    else
        dataset = str2double(input('with-RT filtered EEG dataset to look at: ','s'));
    end 
    
    mbbi_res_flux_est(ALLEEG, EEG, dataset, filepath, dataset_name);

end
