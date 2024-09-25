% jen_logisticregression() - Determine linear discriminating vector between two datasets.
%                            using logistic regression.
%
% THIS IS A MODIFIED VERSION OF ADAM GERSON's FUNCTION
% pop_logisticregression
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Usage:
%   >> jen_logisticregression( ALLEEG, datasetlist, chansubset, chansubset2, trainingwindowlength, trainingwindowoffset, regularize, lambda, lambdasearch, eigvalratio, vinit, show, LOO, domain, PSDofEEG);
%
% Inputs:
%   ALLEEG               - array of datasets
%   datasetlist          - list of datasets
%   chansubset           - vector of channel subset for dataset 1          [1:min(Dset1,Dset2)]
%   chansubset2          - vector of channel subset for dataset 2          [chansubset]
%   trainingwindowlength - Length of training window in samples            [all]
%   trainingwindowoffset - Offset(s) of training window(s) in samples      [1]
%   regularize           - regularize [1|0] -> [yes|no]                    [0]
%   lambda               - regularization constant for weight decay.       [1e-6]
%                            Makes logistic regression into a support 
%                            vector machine for large lambda
%                            (cf. Clay Spence)
%   lambdasearch         - [1|0] -> search for optimal regularization 
%                            constant lambda
%   eigvalratio          - if the data does not fill D-dimensional
%                            space, i.e. rank(x)<D, you should specify 
%                            a minimum eigenvalue ratio relative to the 
%                            largest eigenvalue of the SVD.  All dimensions 
%                            with smaller eigenvalues will be eliminated 
%                            prior to the discrimination. 
%   vinit                - initialization for faster convergence           [zeros(D+1,1)]
%   show                 - display discrimination boundary during training [0]
%   LOO                  - perform leave-one-out cross validation [0]
%   domain               - 'time' if time windows were input
%                          'frequency' if frequency windows were infput
%                           ['time']
%   PSDofEEG             - structure containing power spectral densities of
%                           the 2 datasets and the frequencies at which
%                           they were calculated 
%                               PSDofEEG.PSD 
%                               PSDofEEG.frequency
%                               PSDofEEG.fs
%
% References:
%
% @article{gerson2005,
%       author = {Adam D. Gerson and Lucas C. Parra and Paul Sajda},
%       title = {Cortical Origins of Response Time Variability
%                During Rapid Discrimination of Visual Objects},
%       journal = {{NeuroImage}},
%       year = {in revision}}
%
% @article{parra2005,
%       author = {Lucas C. Parra and Clay D. Spence and Adam Gerson 
%                 and Paul Sajda},
%       title = {Recipes for the Linear Analysis of {EEG}},
%       journal = {{NeuroImage}},
%       year = {in revision}}
%
% Authors: Adam Gerson (adg71@columbia.edu, 2004),
%          with Lucas Parra (parra@ccny.cuny.edu, 2004)
%          and Paul Sajda (ps629@columbia.edu 2004)
%          Modified by Jen Walz (jw2552@columbia.edu 2009)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Adam Gerson, Lucas Parra and Paul Sajda
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Jen NOTE: it may be unnecessary to return the y and a values generated using
% training data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------- %
% Code history   %
% -------------- %

% 5/21/2005, Fixed bug in leave one out script, Adam
% 10/27/2009, Edited for direct call from command window, Jen
% 12/14/2009, Modified to also use power as features, Jen
% 8/20/2010, Commented everything extensively, Jen

% 2013.04.28, added using residual y-value option, by LH.

% Updated by LH, 190108
% Updated residual y-value computation. Used to reshape y-value manually
% via a set number of 50 (samples in each training window), yet apparently
% this should be modified accordingly to reflect samples in windows when
% EEG is recorded w/ different sampling rate

% Note also, that the only difference LR.mat in regressing or not
% regressing out RT, lies in the forward model (a), NOT the y-value itself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ALLEEG, com, Az, Azloomean, weights, bias, y, a, yLOO, aLOO] = jen_logisticregression(ALLEEG, varargin)

%% initializations and checks

bootstrap = 0; % Jen add
Azloomean = []; % JW add so LOO Az can be returned
                
if bootstrap, showaz=0; else showaz=1; end;

com = '';
if nargin < 1
    help jen_logisticregression;
    return;
end;   
if isempty(ALLEEG)
    error('jen_logisticregression(): cannot process empty sets of data');
end;

result = varargin;

%% read input and set defaults

if isempty(result), return; end;

setlist   	 = result{1};
if isempty(setlist), setlist = [1 2]; end

chansubset   = result{2};
if isempty(chansubset), chansubset = 1:min(ALLEEG(setlist(1)).nbchan,ALLEEG(setlist(2)).nbchan); end;

chansubset2  = result{3};
if isempty(chansubset2), chansubset2=chansubset; end

trainingwindowlength    = result{4};
if isempty(trainingwindowlength)
    trainingwindowlength = ALLEEG(setlist(1)).pnts;
end;

trainingwindowoffset    = result{5};
if isempty(trainingwindowoffset)
    trainingwindowoffset = 1;
end;

regularize=result{6};
if isempty(regularize),
    regularize=0;
end

lambda=result{7};
if isempty(lambda),
    lambda=1e-6;
end

lambdasearch=result{8};
if isempty(lambdasearch),
    lambdasearch=0;
end

eigvalratio=result{9};
if isempty(eigvalratio),
    eigvalratio=1e-6;
end

vinit=result{10};
if isempty(vinit),
    vinit=zeros(length(chansubset)+1,1);
end
vinit=vinit(:); % Column vector

show=result{11};
if isempty(show),
    show=0;
end

LOO=result{12};
if isempty(LOO),
    show=0;
end

domain = result{13};
PSDofEEG = result{14};
if isempty(domain) || strcmp(domain,'time'), % we are doing LR in time
    domain = 'time';
elseif strcmp(domain,'frequency'), % we are doing LR in frequency
     if isempty(PSDofEEG)
         error('jen_logisticregression(): you must input the PSD and corresponding frequencies of the data')
     end
     freq = PSDofEEG(1).frequency;
     fs = PSDofEEG(1).fs;
else
    error('jen_logisticregression(): improper domain specified')
end

plotresidual_idx = result{15};  % added by LH, 2013.04.28
rt = result{16};

%% check some stuff

if strcmp(domain,'time'),
    if max(trainingwindowlength+trainingwindowoffset-1)>ALLEEG(setlist(1)).pnts,
        error('jen_logisticregression(): training window exceeds length of dataset 1');
    end
    if max(trainingwindowlength+trainingwindowoffset-1)>ALLEEG(setlist(2)).pnts,
        error('jen_logisticregression(): training window exceeds length of dataset 2');
    end
end
if (length(chansubset)~=length(chansubset2)),
    error('Number of channels from each dataset must be equal.');
end

%% ----------------------- make the labels --------------------------

truth=[zeros(trainingwindowlength.*ALLEEG(setlist(1)).trials,1); ...
    ones(trainingwindowlength.*ALLEEG(setlist(2)).trials,1)];

if bootstrap,
truth=truth(randperm(length(truth)));
end

%% initialize so we can store info in the EEGLAB structure

for i=1:2,
    ALLEEG(setlist(i)).icaweights=zeros(length(trainingwindowoffset),ALLEEG(setlist(i)).nbchan);
    % In case a subset of channels are used, assign unused electrodes in scalp projection to NaN
    ALLEEG(setlist(i)).icawinv=nan.*ones(length(trainingwindowoffset),ALLEEG(setlist(i)).nbchan)';
    ALLEEG(setlist(i)).icasphere=eye(ALLEEG(setlist(i)).nbchan);    
end

%% do it 
    
for i=1:length(trainingwindowoffset),
    
    if strcmp(domain,'time'), % concatenate the datasets
        x=cat(3,ALLEEG(setlist(1)).data(chansubset,...
            trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:), ...
            ALLEEG(setlist(2)).data(chansubset,...
            trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:));
    elseif strcmp(domain,'frequency'), % concatenate the PSDs of the datasets
        x=cat(3,PSDofEEG(1).PSD(chansubset,...
            trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:),...
            PSDofEEG(2).PSD(chansubset,...
            trainingwindowoffset(i):trainingwindowoffset(i)+trainingwindowlength-1,:));
    end
    
    % Rearrange data for logist.m       
	% [D x (T x trials)]' OR [D x (F x trials)]'
    x=x(:,:)';    

	% find the weights by training on all data within that window
    v = logist(x,truth,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
    
    % apply weights to all data within that window
    % this is number of samples per window x number of total trial in both sets
    
    y = x*v(1:end-1) + v(end); 
    
    % probability distribution including all samples
    bp = bernoull(1,y); 
    
    % --------------------------- %
    % Variable dimension          %
    % --------------------------- %
    
    % Below description of variable's dimension added by LH, 190108
    
    % NOTE - this is within the loop for each different window offset, thus
    % the final vector (y) would be of the current dimension * # of windows
    
    % x
    %   | in theory: [samples in training window * (# of std trials + # of odd trials)] * # of channels
    %   | in practice: [25 * (392 + 95)] * 63 | pEf, 180830, Sub35
    
    %   | 50ms training window is equivalent to 25 samples, per a 500Hz sampling rate EEG
    %   | 63 channles since Brain Products' 64 channels include ECG, and 
    %   |             we've rejected that channel prior LR
    
    % v
    %   | in theory: [# of channels] * 1
    %   | in practice: 64 * 1 | pEf, 180830, Sub35
    
    %   | one more row when compared w/ x, for constant setting. The actual
    %   |              rows that multiplies w/ x is of the same dimension as EEG channels, i.e. 63   
    
    % y
    %   | in theory: [samples in training window * (# of std trials + # of odd trials)] * 1
    %   | in practice: [25 * (392 + 95)] * 1 | pEf, 180830, Sub35      

    
    [Az,Ry,Rx] = rocarea(bp,truth); % find Az on training data (using all samples)
    if showaz, % show user the performance for the window
        if strcmp(domain,'time'),
            windowcenter = (trainingwindowoffset(i) + trainingwindowlength/2)/ALLEEG(setlist(1)).srate*1000 + ALLEEG(setlist(1)).xmin*1000;

            fprintf('Window Center: %d ms (offset is %d samples)\n   Az: %6.2f\n', ...
            round(windowcenter), trainingwindowoffset(i), Az); 
        elseif strcmp(domain,'frequency'),
            windowcenter = (trainingwindowoffset(i) + trainingwindowlength/2)/fs + freq(1);

            fprintf('Window Center: %6.2f Hz\n   Az: %6.2f\n', windowcenter, Az); 
        end
    end
    
    % forward model using all samples within window
    
    % added by LH, 2013.04.28
    % updated plotresidual_idx == 2, new if clause by LH, 2013.07.17
        
%     if plotresidual_idx == 1    % regress out RT from both standards and oddballs
%         
%         y_reshape = reshape(y,(ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),trainingwindowlength);
%         residual_y = zeros((ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),trainingwindowlength);
%         
%         for idx = 1:trainingwindowlength
%             
%             p = polyfit(rt,y_reshape(:,idx),1);
%             yfit = polyval(p,rt);
%             
%             residual_y(:,idx) = y_reshape(:,idx)-yfit;
%             
%         end
%         
%         a = reshape(residual_y,trainingwindowlength*(ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),1) \ x;
%         
%     elseif plotresidual_idx == 2    % ONLY regress out RT from oddballs

    % --------------------------- %
    % Variable dimension          %
    % --------------------------- %
    
    % Below description of variable's dimension added by LH, 190108  

    % NOTE - this is within the loop for each different window offset, thus
    % the final vector (y) would be of the current dimension * # of windows    
    
	% y
    %   | in theory: [samples in training window * (# of std trials + # of odd trials)] * 1
    %   | in practice: [25 * (392 + 95)] * 1 | pEf, 180830, Sub35   

	% y_reshape
    %   | in theory: [(# of std trials + # of odd trials)] * samples in training window
    %   | in practice: [(392 + 95)] * 25 | pEf, 180830, Sub35  
    
    % for each sample in each window
    %                | y_reshape(odd_trial,sample): [# of odd trials] * 1
    %                | rt:                          [# of odd trials] * 1
    
    %                | fit rt to y_reshape(odd_trial,sample), get p, then yfit
    
    %                | residual_y(std_trial,sample): y_reshape(std_trial,sample)
    %                | residual_y(odd_trial,sample): y_reshape(odd_trial,sample) - yfit
    
	% a
    %   | in theory: [# of channels] * 1
    %   | in practice: 63 * 1 | pEf, 180830, Sub35  
    
        
    %%% Following reshape() expression updated by LH, 190108 %%%
    
    % Use trainingwindowlength to reshape, instead of a fixed number of 50
    % (okay in Plos One paper as we have 1kHz sampling rate of EEG, whilst
    % it'll be an issue for pEf analysis as we have 500Hz sr)

    if plotresidual_idx == 2    % ONLY regress out RT from oddballs
        
        %%% DEPRECATED USAGE 
        % y_reshape = reshape(y,(ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),50);
        
        % y_reshape: (# of std trials + # of odd trials) * samples in each training window        
        
        %%% In-use        
        
        y_reshape = reshape(y,(ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),trainingwindowlength);
        residual_y = zeros((ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),trainingwindowlength);
        
        for idx = 1:trainingwindowlength
            
            p = polyfit(rt,y_reshape((ALLEEG(setlist(1)).trials+1):end,idx),1);
            yfit = polyval(p,rt);
            
            residual_y(1:ALLEEG(setlist(1)).trials,idx) = y_reshape(1:ALLEEG(setlist(1)).trials,idx);   % standard y's remain the same
            residual_y((ALLEEG(setlist(1)).trials+1):end,idx) = y_reshape((ALLEEG(setlist(1)).trials+1):end,idx)-yfit;  % regress out RT from oddball y's
            
        end
        
        a = reshape(residual_y,trainingwindowlength*(ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials),1) \ x;        
        
    elseif plotresidual_idx == 0
        
        a = y \ x;
    
    end
    
    % JW add:
    % find only one y per trial (to return)
    ymean(:,i) = mean(reshape(y, ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials, trainingwindowlength),2);
    % save the a's to return
    a2return(:,i) = a'; 
    % save the weights to return
    weights(:,i) = v(1:end-1)';
    bias(i) = v(end); 
    
    % old stuff from Adam:
    
%     asetlist1=y(1:trainingwindowlength.*ALLEEG(setlist(1)).trials) \ x(1:trainingwindowlength.*ALLEEG(setlist(1)).trials,:);
%     asetlist2=y(trainingwindowlength.*ALLEEG(setlist(1)).trials+1:end) \ x(trainingwindowlength.*ALLEEG(setlist(1)).trials+1:end,:);

    % save the weights as ICA weights in ALLEEG structure 
    ALLEEG(setlist(1)).icaweights(i,chansubset)=v(1:end-1)';
    ALLEEG(setlist(2)).icaweights(i,chansubset2)=v(1:end-1)';
    
    ALLEEG(setlist(1)).icawinv(chansubset,i)=a'; % consider replacing with asetlist1
    ALLEEG(setlist(2)).icawinv(chansubset2,i)=a'; % consider replacing with asetlist2
    
    if LOO
        
        clear y;
        
        N=ALLEEG(setlist(1)).trials+ALLEEG(setlist(2)).trials; % total # trials
        %ploo=zeros(N*trainingwindowlength,1); % initialize LR probability

        parfor looi=1:length(x)./trainingwindowlength, % for each trial
            
            indx=ones(N*trainingwindowlength,1);    % ones vector length of data
            
            % select indices corresponding to trial and set to 0
            indx((looi-1)*trainingwindowlength+1:looi*trainingwindowlength)=0;
            tmp = x(find(indx),:);          % LOO data
            tmpt = [truth(find(indx))];     % truth labels for LOO data
            
            % get weights using all LOO data
            vloo(:,looi)=logist(tmp,tmpt,vinit,show,regularize,lambda,lambdasearch,eigvalratio);
            
            % apply to the left out data to yet y's for each sample in that trial
            y(:,looi) = [x((looi-1)*trainingwindowlength+1:looi*trainingwindowlength,:) ones(trainingwindowlength,1)]*vloo(:,looi);
            
            % average y's to obtain one y for the trial
            yloomean(looi,i)=mean(y(:,looi));
            
            % find probability distribution using y's for all samples in that trial
            %ploo((looi-1)*trainingwindowlength+1:looi*trainingwindowlength) = bernoull(1,y(:,looi));
            % find probability of trial belonging to class 
            ploomean(looi)=bernoull(1,yloomean(looi,i));            
        end
        
        % truth labels - just one for each trial
        truthmean=([zeros(ALLEEG(setlist(1)).trials,1); ones(ALLEEG(setlist(2)).trials,1)]);
        
        % find Az from probability distribution - one value per sample
        %[Azloo,Ryloo,Rxloo] = rocarea(ploo,truth);
        % find Az from probability distribution - one value per trial
        [Azloomean,Ryloomean,Rxloomean] = rocarea(ploomean,truthmean);
        fprintf('   LOO Az: %6.2f\n',Azloomean);

        % JW add - forward model generated using LOO y values for all samples

        %%% Following reshape() expression updated by LH, 190108 %%%

        % Use trainingwindowlength to reshape, instead of a fixed number of 50
        % (okay in Plos One paper as we have 1kHz sampling rate of EEG, whilst
        % it'll be an issue for pEf analysis as we have 500Hz sr)        
        
%         if plotresidual_idx == 1    % regress out RT from standards+oddballs
% 
%             residual_y = zeros(N,trainingwindowlength);
%             y_temp = y';
% 
%             for idx = 1:size(y_temp,2)
% 
%                 p = polyfit(rt,y_temp(:,idx),1);
%                 yfit = polyval(p,rt);
% 
%                 residual_y(:,idx) = y_temp(:,idx)-yfit;
% 
%             end
% 
%             a = reshape(y_temp,size(y_temp,1)*size(y_temp,2),1) \ x;
%             aLOO(:,i) = a';
% 
%         elseif plotresidual_idx == 2    % regress out RT from oddballs

        % --------------------------- %
        % Variable dimension          %
        % --------------------------- %

        % Below description of variable's dimension added by LH, 190108  

        % NOTE - this is within the loop for each different window offset, thus
        % the final vector (y) would be of the current dimension * # of windows    

        % y
        %   | in theory: [samples in training window] * (# of std trials + # of odd trials)
        %   | in practice: [25] * (392 + 95) | pEf, 180830, Sub35   

        % y_temp, i.e. y'
        %   | in theory: [(# of std trials + # of odd trials)] * samples in training window
        %   | in practice: [(392 + 95)] * 25 | pEf, 180830, Sub35   
        
        
        % for each sample in each window
        %                | y_temp(odd_trial,sample): [# of odd trials] * 1
        %                | rt:                       [# of odd trials] * 1

        %                | fit rt to y_temp(odd_trial,sample), get p, then yfit

        %                | residual_y(std_trial,sample): y_temp(std_trial,sample)
        %                | residual_y(odd_trial,sample): y_temp(odd_trial,sample) - yfit

        % a
        %   | in theory: [# of channels] * 1
        %   | in practice: 63 * 1 | pEf, 180830, Sub35 

        if plotresidual_idx == 2    % regress out RT from oddballs
        
            residual_y = zeros(N,trainingwindowlength);
            y_temp = y';

            for idx = 1:trainingwindowlength

                p = polyfit(rt,y_temp((ALLEEG(setlist(1)).trials+1):end,idx),1);
                yfit = polyval(p,rt);

                residual_y(1:ALLEEG(setlist(1)).trials,idx) = y_temp(1:ALLEEG(setlist(1)).trials,idx);
                residual_y((ALLEEG(setlist(1)).trials+1):end,idx) = y_temp((ALLEEG(setlist(1)).trials+1):end,idx)-yfit;

            end

            a = reshape(y_temp,size(y_temp,1)*size(y_temp,2),1) \ x;
            aLOO(:,i) = a';
            
        elseif plotresidual_idx == 0

            a = reshape(y,size(y,1)*size(y,2),1) \ x;
            aLOO(:,i) = a'; 

        end
        
    end
    
end 

   
clear y
% JW add:
% store y's for each set separately
y{1} = ymean(1:ALLEEG(setlist(1)).trials,:); 
y{2} = ymean(ALLEEG(setlist(1)).trials+1:end,:);
if LOO,
yLOO{1} = yloomean(1:ALLEEG(setlist(1)).trials,:); 
yLOO{2} = yloomean(ALLEEG(setlist(1)).trials+1:end,:);
end
a = a2return;

%% -- recompute activations (temporal sources) and check dataset integrity

eeg_options;
for i=1:2, % for each dataset
    if option_computeica
        if strcmp(domain,'time'),
            ALLEEG(setlist(i)).icaact    = (ALLEEG(setlist(i)).icaweights*ALLEEG(setlist(i)).icasphere)*reshape(ALLEEG(setlist(i)).data, ALLEEG(setlist(i)).nbchan, ALLEEG(setlist(i)).trials*ALLEEG(setlist(i)).pnts);
            ALLEEG(setlist(i)).icaact    = reshape( ALLEEG(setlist(i)).icaact, size(ALLEEG(setlist(i)).icaact,1), ALLEEG(setlist(i)).pnts, ALLEEG(setlist(i)).trials);
        elseif strcmp(domain,'frequency'),
            ALLEEG(setlist(i)).icaact    = (ALLEEG(setlist(i)).icaweights*ALLEEG(setlist(i)).icasphere)*reshape(PSDofEEG(i).PSD, ALLEEG(setlist(i)).nbchan, ALLEEG(setlist(i)).trials*length(freq));
            ALLEEG(setlist(i)).icaact    = reshape( ALLEEG(setlist(i)).icaact, size(ALLEEG(setlist(i)).icaact,1), length(freq), ALLEEG(setlist(i)).trials);
        end
    end    
end

%% save the datasets

% this will store the LR weights as ICA weights so scalp maps can be
% plotted later
for i=1:2, [ALLEEG]=eeg_store(ALLEEG,ALLEEG(setlist(i)),setlist(i)); end

com = sprintf('jen_logisticregression( %s, [%s], [%s], [%s], [%s]);', inputname(1), num2str(setlist), num2str(chansubset), num2str(trainingwindowlength), num2str(trainingwindowoffset));


%fprintf('\nLR done.\n');

