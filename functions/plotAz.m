%
% [maxAz,maxwindow,linehandle] = plotAz(LR,plotcolor,fighandle,filepath,savefig)
%
% Plot and save Az vs windowcenter (time or frequency window)
% Save can be turned off if you are using this function to overlay plots
% and do not want figure saved after each overlay.
%
% INPUT:
% LR = structure containing LR results (see runLR())
%
% OPTIONAL INPUT:
% plotcolor     = color of plot, useful if overlaying conditions
%                 cell of strings {'Training_plotcolor'; 'LOO_plotcolor'}
%                 or if single string entered, plotcolor will be used for both
%                 AzTraining and AzLOO
%                 [] for defaults: 'b' for training and 'g' for LOO
% fighandle     = handle to figure on which you want to plot
%                 [] for new figure
% filepath      = directory in which to save figures, default is cd
% savefig       = [0 | 1] default 1 (save the figure)
%
% OUTPUT:
% maxAz         = max LOO Az value (max mean LOO Az value if AzLOO is matrix)
% maxwindow     = center of window with max Az
% linehandle    = handle to the line plotted, useful if overlaying conditions
%                 if both AzTraining and AzLOO are plotted, handle returned
%                 will be LOO line
%
% Written by Jen Walz 8/2/10 (jw2552@columbia.edu)

function [maxAz,maxwindow,linehandle] = plotAz(LR,plotcolor,fighandle,filepath,savefig)

%% defaults

% (LR, varargin)
% for i=1:2:length(varargin)
    
if nargin<5 || isempty(savefig), savefig = 1; end
if nargin<4 || isempty(filepath), filepath = cd; end
% make new fig and get handle for it
if nargin<3 || isempty(fighandle), fighandle = figure; end
% define plot colors
if nargin<2 || isempty(plotcolor), 
    plotcolor = {'b'; [0 0.5 0]}; % {training; LOO}
elseif ischar(plotcolor),
    plotcolor = {plotcolor; plotcolor};
end
maxAz = [];
maxwindow = [];


%% initializations

if strcmp(LR.domain,'time'), unit = 'ms'; else unit = 'Hz'; end
%exptname = [LR.epochtype ' ' LR.subjectID ' ' LR.condition ' ' num2str(LR.windowlength) unit 'Windows'];

%% find mean across subjects if required

%if size(LR.AzTraining,1)==1 && size(LR.AzLOO,1)==1, 
    %filepath = [filepath '/Az_plots']; 
%end
if size(LR.AzLOO,1)>1, % if a matrix was input
    % find mean LOO Az values across subjects for all windows
    stderrLOO = std(LR.AzLOO,1)/sqrt(size(LR.AzLOO,1)); % also find the std error
    LR.AzLOO = mean(LR.AzLOO,1);   
    %filepath = [filepath '/Mean_Az_Plots'];
end
if size(LR.AzTraining,1)>1,
    stderrTraining = std(LR.AzTraining,1)/sqrt(size(LR.AzTraining,1)); % also find the std error
    LR.AzTraining = mean(LR.AzTraining,1);
    %filepath = [filepath '/Mean_Az_Plots'];
end
if ~exist(filepath,'dir'), mkdir(filepath), end

%% display max LOO Az value to user

if ~isempty(LR.AzLOO),
% find max LOO Az and the index of window with max Az LOO
winidx = find(LR.AzLOO==max(LR.AzLOO));
maxAz = LR.AzLOO(winidx);   % find max mean Az
maxwindow = LR.windowcenter(winidx); % find that window center

if exist('stderrLOO','var'),
    fprintf('\nThe maximum mean LOO Az is %6.2f with std error %6.2f at %6.2f %s\n', ...
        maxAz,stderrLOO(winidx),maxwindow,unit)
else
    fprintf('\nThe maximum LOO Az is %6.2f at %6.2f %s\n', ...
        maxAz,maxwindow,unit)   
end
end

%% generate figure

fprintf(['\nPlotting LR Classifier Performance... '])
       
figure(fighandle)
hold on

if ~isempty(LR.AzTraining),
    linehandle = plot(LR.windowcenter,LR.AzTraining,'Color',plotcolor{1});
end
if ~isempty(LR.AzLOO),
    % plot LOO Az and return handle for 
    linehandle = plot(LR.windowcenter,LR.AzLOO,'Color',plotcolor{2});
end

if ~isempty(LR.AzTraining) && ~isempty(LR.AzLOO),
    % add legend if we are plotting both training and LOO
    legend('Training','LOO')
end

if ~isempty(LR.AzTraining) && exist('stderrTraining','var'),
    % if we are plotting a subject mean, add standard error patches
    fill([LR.windowcenter fliplr(LR.windowcenter)], ...
        [LR.AzTraining+stderrTraining fliplr(LR.AzTraining-stderrTraining)],plotcolor{1},...
        'FaceAlpha',0.2,'EdgeColor','none')
end
if ~isempty(LR.AzLOO) && exist('stderrLOO','var'),
    % if we are plotting a subject mean, add standard error patches
    fill([LR.windowcenter fliplr(LR.windowcenter)], ...
        [LR.AzLOO+stderrLOO fliplr(LR.AzLOO-stderrLOO)],plotcolor{2},...
        'FaceAlpha',0.2,'EdgeColor','none')
end

plot(LR.windowcenter,0.5*ones(length(LR.windowcenter)),'k--') % show chance line
plot(LR.windowcenter,0.75*ones(length(LR.windowcenter)),'k:') % show significance line
plot(zeros(1,2),[0.3 1],'Color',[0.8 0.8 0.8]) % draw line at zero
axis([min(LR.windowcenter)-LR.windowlength/2 max(LR.windowcenter)+LR.windowlength/2 0.3 1])

% label figure

xlabel(['Center of Window [' unit ']'],'FontName','Geneva');%,'FontSize',14)
title(['LR Classifier Performance for ' LR.subjectID(LR.subjectID~='_') ' ' ...
    LR.epochtype ' ' LR.sessiontype ' ' LR.condition ' using ' ...
    num2str(LR.windowlength) ' ' unit ' Windows'],...
    'FontName','Geneva');%,'FontSize',14)
ylabel('Area Under ROC Curve','FontName','Geneva');%,'FontSize',14)
set(gcf, 'Position', [0 0 1200 800]);     % maximize figure

%% save figure

if savefig,
figname = [filepath '/Azplot_' LR.subjectID '_' ...
    LR.epochtype '_' LR.sessiontype '_' LR.condition];
figname(figname==' ')='_';
    
% as a safety, don't overwrite existing file. put in new folder instead
if exist([filepath '/' figname '.fig'],'file'),
    time = datestr(clock); time(time==':')='-'; time(time==' ')='_';
    newpath = [filepath '/moved_here_' time];
    mkdir(newpath)
    movefile([filepath '/Azplot_' LR.subjectID '*'],newpath);
end

% save the figure to the current directory with name corresponding to the
% saved result;
saveas(gcf, figname, 'fig'),
saveprettyfig(figname),
end

fprintf(['Done.\n'])
