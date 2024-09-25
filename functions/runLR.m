% runLR() runs the logistic regression plugin using jen_logisticregression()
%
% USAGE:
%
% [ALLEEG, LR, LRfilename] = runLR(ALLEEG,setlist,chansubset,subjectID,condition, ...
%                           domain,windowlength,windowcenter,LOO,filepath,sessiontype);
%
% 
% INPUT:
%
% ALLEEG       = EEG data structure from EEGLAB
% setlist      = numbers of the EEGLAB datasets to classify
%                   default [1 2]
%                   (note order will be reversed before input to LR function)
% chansubset   = channels to use in classification, default all []
% subjectID    = string identifying the subject
% condition    = string specifying the discrimination condition
% domain       = string determining the domain in which we will classify
%                ['time' | 'frequency'] 
%                default is 'time'
% windowlength = length of training window in ms or Hz depending on domain
%                   default 50 ms or 4 Hz
% windowcenter = scalar or vector of training window centers in ms or Hz
%                locked to time of event or zero frequency 
%                   default 0 
% LOO          = on (1) or off (0), default off
% filepath     = string name of directory in which to save results
%                defaults to current directory
% sessiontype  = if there are multiple types of runs within the same
%                   experiment, use this for titles 
%                   (ex: visual or auditory oddball sessions)
%
% OUTPUT:
%
% ALLEEG       = EEG data structure
% LR           = struct containing results of logistic regression
% LRfilename   = name of file in which results are saved (saved in filepath)
%
% This function calls jen_logisticregression.m, which was modified
% from pop_logisticregression.m (by Adam Gerson) to 
% 1. run from command line
% 2. have option to classify in frequency domain.
% 3. save and return results
% 4. print center of window in ms or Hz instead of onset in samples
% 
% see also: jen_logisticregression()
%

% created by Jen Walz, 12/14/09
% last update JMW 10/22/10

% 2013.04.28, added using residual y-value option, by LH.


function [ALLEEG, LR, LRfilename] = runLR(varargin)

%% --------------------- initial check --------------------------

ALLEEG = varargin{1};
    
LR = struct;

%% -------------------- set defaults ---------------------------

if nargin>1, setlist = varargin{2}; end
if ~exist('setlist','var')||isempty(setlist), setlist = [1 2]; end       

if nargin>2, LR.chansubset = varargin{3}; end
if ~isfield(LR,'chansubset')||isempty(LR.chansubset), LR.chansubset = []; end       

if nargin>3, LR.subjectID = varargin{4}; end
if ~isfield(LR,'subjectID')||isempty(LR.subjectID), LR.subjectID = 'One Subject'; end 

if nargin>4, LR.condition = varargin{5}; end
if ~isfield(LR,'condition')||isempty(LR.condition), LR.condition = ''; 
else
condition4fname = ['_' LR.condition];
condition4fname(condition4fname==' ')='_';
end

if nargin>5, LR.domain = varargin{6}; end
if ~isfield(LR,'domain')||isempty(LR.domain), LR.domain = 'time'; end     

if nargin>6, LR.windowlength = varargin{7}; end
if nargin>7, LR.windowcenter = varargin{8}; end

if nargin>8, LOO = varargin{9}; end
if ~exist('LOO','var') || isempty(LOO), LOO = 1; end

if nargin>9, filepath = varargin{10}; end
if ~exist('filepath','var') || isempty(filepath), filepath = cd; end

if nargin>10, LR.epochtype = varargin{11}; end
if ~isfield(LR,'epochtype') || isempty(LR.epochtype),
    LR.epochtype = 'Stimulus Locked';
end
epochtype4fname = ['_' LR.epochtype];
epochtype4fname(epochtype4fname==' ')='';

if nargin>11, LR.sessiontype = varargin{12}; end
if ~isfield(LR,'sessiontype') || isempty(LR.sessiontype), 
    LR.sessiontype = ''; 
end

if nargin>12, LR.plotresidual = varargin{13}; end
if ~isfield(LR,'plotresidual') || isempty(LR.plotresidual), 
    LR.plotresidual = 0; 
end

if nargin>13, LR.rt = varargin{14}; end
if ~isfield(LR,'rt') || isempty(LR.rt),
    LR.rt = [];
end
% added by LH, 2013.04.28

% initializations dependent upon domain choice
switch LR.domain
    case 'time'
        if ~isfield(LR,'windowlength')||isempty(LR.windowlength), LR.windowlength = 50; end    
        if ~isfield(LR,'windowcenter')||isempty(LR.windowcenter), LR.windowcenter = 0; end    
        
        fs = ALLEEG(setlist(1)).srate/1000;     % samples per ms
        shift = ALLEEG(setlist(1)).xmin*1000;   % epoch start relative to stim in ms
        totallength = ALLEEG(setlist(1)).pnts;  % samples per epoch
        PSDofEEG = [];
        unit = 'ms';
        
        fprintf(['\nPerforming Time Domain Analysis\n'])
        
    case 'frequency'
        if ~isfield(LR,'windowlength')||isempty(LR.windowlength), LR.windowlength = 4; end    
        if ~isfield(LR,'windowcenter')||isempty(LR.windowcenter), LR.windowcenter = 0; end    
        
        fs = 5;  % samples per Hz               
        % specify frequency range for PSD calculation
        frequency = 0:1/fs:50; % in Hz - This was initially up to 100 Hz
        shift = frequency(1); % start of PSD relative to 0 freq in Hz
        totallength = length(frequency); % samples in PSD
        PSDofEEG = struct;
        unit = 'Hz';
        
        % directory for saving PSDs
        PSDfilepath = [filepath '/PSD_data'];
        if ~exist(PSDfilepath,'dir'), mkdir(filepath,'PSD_data'), end
        method = 'Welch';
        fprintf(['\nPerforming Frequency Domain Analysis using ' method ' Method...\n'])
end


%% -------- input to LR plugin function jen_logisticregression() --------

% standard inputs
%chansubset = [];
%chansubset2 = [];
regularize = 1;     % might as well try regularization
lambda = 1e-6;      % try 1e-4 and 1e-3 ?
lambdasearch = 0;   % search for the optimal regularization constant
eigvalratio = 1e-6; % 1e-6 is default, only increase if rank deficient
vinit = [];
show = 0;           % don't show the discrimination boundary

%% -------------- PSD initializations -------------


if strcmp(LR.domain,'frequency'),
    
    for i=1:2, % for each of the datasets we are classifying
        % first check for existence of PSD 
        fname = ['PSD_' method '_' ALLEEG(setlist(i)).setname];
        PSD = struct;
        if exist([PSDfilepath '/' fname '.mat'],'file'),    
            
            % load saved PSD and add to structure
            % will contain 3 variables: frequency, fs, and PSD
            load([PSDfilepath '/' fname '.mat']);
            
            % check consistency
            if PSD.fs~=fs || PSD.frequency(1)~=frequency(1) || PSD.frequency(end)<frequency(end),
                fprintf(['\nThere is an inconsistency with the current frequency vector\n' ...
                    'and the saved PSD frequencies for %s. \n'...
                'PSD will be recalculated.\n'],ALLEEG(setlist(i)).setname)
                % rename this PSD so we can re-calculate and save
                save([PSDfilepath '/' fname '_' num2str(PSD.fs) 'sampleperHz.mat'],'PSD');
                
                PSD.PSD = []; % to force re-calculation of PSD
            else
            % if we already calculated and saved the PSDs in the cd don't do it again
            fprintf('\nUsing previously saved PSD for DATASET %s.', ALLEEG(setlist(i)).setname)
            end
        end
        if ~isfield(PSD,'PSD') || isempty(PSD.PSD),
            % if we don't already have the PSD stored we need to
            % calculate (or re-calculate) it
            fprintf(['\nCalculating PSD for %s...\n'], ALLEEG(setlist(i)).setname)
            
            [PSD] = createPSD4LR(ALLEEG(setlist(i)),frequency,PSDfilepath,method);
            
        end 
        % create PSDofEEG as struct - this is input to LR plugin
            PSDofEEG(i).PSD = PSD.PSD;
            PSDofEEG(i).frequency = PSD.frequency;
            PSDofEEG(i).fs = PSD.fs;
       
        clear PSD
    end   
    
    % NORMALIZE the PSD
%     
%     if normalizePSD,
%         fprintf('\n\nNormalizing the PSD...')
%         for i=1:2, % each dataset
%             for tr=1:size(PSDofEEG(i).PSD, 3),
%                 for ch=1:size(PSDofEEG(i).PSD, 1),
%                     PSDofEEG(i).PSD(ch,:,tr) = PSDofEEG(i).PSD(ch,:,tr) .* PSDofEEG(i).frequency;
%                 end
%             end
%         end
%         fprintf(' Done.\n')
%     end
end


%% ------------ ADJUST INPUT for LR ----------------------

% find onset of windows from list of centers of windows
% in time domain this is relative to event lock (stim or resp)
windowstart = LR.windowcenter - LR.windowlength/2; % ms or Hz
% adjust window so it corresponds to indices
LR.windowlength_idx = round(LR.windowlength*fs);
LR.windowstart_idx = round((windowstart-shift)*fs);

% ensure we don't exceed limits of epoch or PSD
goodidx = find(LR.windowstart_idx <= totallength-LR.windowlength_idx ...
    & LR.windowstart_idx > 0);
LR.windowstart_idx = LR.windowstart_idx(goodidx);
LR.windowcenter = LR.windowcenter(goodidx);


%% ------------- RUN LOGISTIC REGRESSION -------------------

fprintf(['\nRunning Logistic Regression...\n'])

LR.y{1}=[]; LR.y{2}=[]; LR.yLOO{1}=[]; LR.yLOO{2}=[];
clear y yLOO
for i = 1:length(LR.windowstart_idx), % for each offset
    % do it
    [ALLEEG,com,LR.AzTraining(i),LR.AzLOO(i),LR.weights{i},LR.bias(i), ...
        y{i},LR.a{i},yLOO{i},LR.aLOO{i}] = ...
        jen_logisticregression(ALLEEG, fliplr(setlist), LR.chansubset, LR.chansubset, ...
        LR.windowlength_idx, LR.windowstart_idx(i), regularize, lambda, ...
        lambdasearch, eigvalratio, vinit, show, LOO, LR.domain, PSDofEEG, LR.plotresidual, LR.rt); % added LR.plotresidual option by LH.
    
    % y and yLOO values stored in 2-element cell arrays
    % where y{1} is a matrix size [NtrialsSet1 x Nwindows]
    % and y{2} is a matrix size [NtrialsSet2 x Nwindows]
    LR.y{1} = [LR.y{1} y{i}{2}]; 
    LR.y{2} = [LR.y{2} y{i}{1}];
    LR.yLOO{1} = [LR.yLOO{1} yLOO{i}{2}]; 
    LR.yLOO{2} = [LR.yLOO{2} yLOO{i}{1}];
end

% matrices are easier to work with (size is Nchansubset x Nwindows)
LR.aLOO = cell2mat(LR.aLOO); 
LR.a = cell2mat(LR.a); 
LR.weights = cell2mat(LR.weights);
LR.info = sprintf(['weights, bias, y, and a are from training data:\n'...
            '   y = weights * data + bias  (where y is the mean y for each trial)\n'...
            '   a = y \\ data  (where y has length #Trials*#samplesPerTrial)\n'...
            'yLOO and aLOO are generated using LOO cross validation:\n'...
            '   yLOO = weightsLOO * data + biasLOO  (where yLOO is mean yLOO across all samples in trial)\n'...
            '   (weightsLOO and biasLOO are different for each trial)\n'...
            '   aLOO = yLOO \\ data  (where yLOO has length #Trials*#samplesPerTrial)\n']);% '...
%             'In the case of EEG/fMRI data:\n'...
%             '   If present, aLOO_b4reref and weights_b4reref are as described above, and \n'...
%             '   aLOO is the forward model re-referenced from bipolar pairs\n'...
%             '   to the electrode space (aLOO=A*aLOO_b4reref).\n'...
%             '   weights includes zeros for bipolar pair channels excluded from classification.\n']);
LR.datasetnames{1} = ALLEEG(setlist(1)).filename;
LR.datasetnames{2} = ALLEEG(setlist(2)).filename;
LR.chanlocs = ALLEEG(setlist(1)).chanlocs;

% For 111026 Lisa' experiment, not really have inserted responses for
% standard stimuli, hence leave them right there.

% % store RT in ms - note this RT is based on an assumption - see guessRTs()
% LR.RT{1} = guessRTs(ALLEEG(setlist(1)));
% LR.RT{2} = guessRTs(ALLEEG(setlist(2)));

% show user the max Az values we found
if length(LR.windowstart_idx)>1,
fprintf('\nMax Az = %1.4f at %4.1f %s (%3.1f %s window)\n', ...
    max(LR.AzTraining), LR.windowcenter(find(LR.AzTraining==max(LR.AzTraining),1)), unit, LR.windowlength, unit)
if LOO,
    fprintf('Max LOO Az = %1.4f at %4.1f %s (%3.1f %s window)\n', ...
        max(LR.AzLOO), LR.windowcenter(find(LR.AzLOO==max(LR.AzLOO),1)), unit, LR.windowlength, unit)
end
end

%% --------------------------- SAVE RESULTS ------------------------------

fprintf(['\nSaving Results of Logistic Regression\n'])

LRfilepath = [filepath '/LR_results'];
if ~exist(LRfilepath,'dir'), mkdir(filepath,'LR_results'); end

% save the results overwriting any old results
LRfilename = ['LRresults_' LR.subjectID '_' LR.sessiontype condition4fname ...
    epochtype4fname '.mat'];

% as a safety, don't overwrite existing file. put in new folder instead
if exist([LRfilepath '/' LRfilename],'file'),
    time = datestr(clock); time(time==':')='-'; time(time==' ')='_';
    newpath = [LRfilepath '/moved_here_' time];
    mkdir(newpath)
    movefile([LRfilepath '/' 'LRresults_' LR.subjectID '*'],newpath);
end

save([LRfilepath '/' LRfilename], 'LR')

fprintf(['\nLogistic Regression DONE\n'])

% EOF
