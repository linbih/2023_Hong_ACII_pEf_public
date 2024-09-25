%
% scalpplotseries(LR,windows2map,filepath)
%
% plots and saves a figure showing forward model of discriminating
% component for a series of logistic regression training windows.  
% Forward model that is plotted is from the LOO cross validation.
%
% INPUT
%
% LR          = structure containing results of LR classification
% windows2map    = vector of centers of training windows which we want 
%               to include in map series 
% filepath    = string name of directory in which to save figs,
%               defaults to cd
%
% A scalp map series is generated and saved as a figure 
% in a 'Scalp_Plots' folder within the current directory.
%

% Written by JMW 8/2/10 (jw2552@columbia.edu)
% Last Updata:  9/20/10 JMW

function scalpplotseries(LR,windows2map,filepath)


%% initializations

if nargin<3 || isempty(filepath),
    filepath = cd;
end

% for labels
if strcmp(LR.domain,'time'), unit = 'ms'; else unit = 'Hz'; end

% a is channels x windows for a single subject, or subject x channels x
% windows if we are plotting mean across subjects
a = LR.aLOO; 
%a = LR.a; 

% if iscell(a), a = cell2mat(a); end

if length(size(a))==3, % a is size channels x windows x subjects
    % we want to plot a series of mean projections across subs
    
    filepath = [filepath '/Mean_Scalp_Plots'];
    if ~exist(filepath,'dir'), mkdir(filepath); end

    collimits = [-0.5 0.5]; % because we normalize the forward models
    %collimits = [-4 4]; % if we don't normalize the forward models
    
    % normalize each subject 
    for sub = 1:size(a,3),
        a(:,:,sub) = a(:,:,sub) / max(max(abs(a(:,:,sub))));
    end
    
    % find the mean
    a = squeeze(mean(a,3));
    
else
    filepath = [filepath '/Scalp_Plots'];
    if ~exist(filepath,'dir'), mkdir(filepath); end
    
    collimits = [-3 3]; %used to be -5 5
end

% topoplot options
typeplot = 0;   % plot component maps
items = 1;      % for topoplot function call
rowscols = [];  % near square
plotdip = 0;    % no associated dipoles
chans2plot = [];  % plot all chans
colmap = [];    % use default colormap

% we can only plot scalp maps at latencies for which we have weights
windows2map = intersect(windows2map,LR.windowcenter);

% for display: scalp plot figure no more than 9 heads across
ncols = min(length(windows2map),9); nrows = ceil(length(windows2map)/9);

% check to make sure the windows selected for maps are the windows
% for which we have logistic regression results
if isempty(windows2map),
    fprintf(['\nThere is no data available for the windows ' ...
        'which you selected for scalp maps.' ...
        '\nNo scalp projections will be plotted.\n'])
    return
end


%% do it! plot some scalp projections  :-)

fprintf(['\nPlotting (%s domain) Component Scalp Projections ' ...
    'for %s %s...\n'], ...
    LR.domain, LR.subjectID, LR.condition)

% initialize and name figure
figure('Color',[1 1 1],'NumberTitle','off', ...
    'Name',['Discriminant Component Scalp Projections - ' ...
    LR.subjectID ' - ' LR.condition], ...
    'Position',[100 100 1400 800])

% for each of the windows we selected for component maps
for mapidx = 1:length(windows2map),
    % get the index of the window we are plotting
    winidx = find(LR.windowcenter==windows2map(mapidx));
    h = subplot(nrows,ncols,mapidx);
    
    % plot one scalp map
    topoplot( a(:,winidx), LR.chanlocs, ...
        'headrad', 0.5, 'maplimits', collimits, 'style', 'map', 'electrodes', 'off', ...
        'plotchans', chans2plot, 'verbose', 'off', 'drawaxis', 'off', ...
        'hcolor', 'k', 'whitebk', 'on');

    % label the subplot with the window center and Az
    title({[num2str(windows2map(mapidx)) ' ' unit];['Az = ' num2str(LR.AzLOO(winidx),2)]},...
        'FontName','Geneva')
    
    % total workaround to make a big title on the bottom
    % JEN: fix this (look at make_mean_figs_like_paper)
    if mapidx == (nrows-1)*ncols+ceil(ncols/2),
        text('String',['Discriminant Component Scalp Projections - ' ...
     LR.subjectID(LR.subjectID~='_') ' - ' LR.sessiontype ' ' LR.condition ' - ' LR.epochtype],...
     'Units','normalized','VerticalAlignment','middle',...
     'HorizontalAlignment','center',...
     'Position',[0.5, -0.5],'FontSize',14,'FontName','Geneva','Color','k');
    end
    
end

% plot a colorbar
pos = get(h,'Position');
cbar = colorbar;
%set(cbar,'Position',[pos(1)+pos(3)+0.03, pos(2)+0.04, 0.01, pos(4)-0.07]);
set(cbar,'Position',[pos(1)+pos(3)+0.04, 0.15, 0.01, 0.75]);

set(gcf,'Color',[1 1 1])    
 
%  save the figure

% file names for scalp map figures
fname = ['ScalpMaps_' LR.subjectID '_' LR.sessiontype '_' LR.condition '_' LR.epochtype];
fname = fname(fname~=' ');
% as a safety, don't overwrite existing file. put in new folder instead
if exist([filepath '/' fname '.fig'],'file') && exist([filepath '/' fname '.jpg'],'file'),
    now = datestr(clock); now(now==':')='-';
    newpath = [filepath '/moved_here_' now];
    mkdir(newpath)
    movefile([filepath '/ScalpMaps_' LR.subjectID '*'],newpath);
end

% save the scalp map figure
figname = [filepath '/' fname];
saveas(gcf, figname, 'fig'),
saveprettyfig(figname),

