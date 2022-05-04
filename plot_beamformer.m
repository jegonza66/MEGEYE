clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
addpath('markus/')
experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319'};

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

cfg.myfname = ['dataEpoched_',sessionfilenames{isubject},'_src413nofilt_merged.mat'];
load([cfg.directory,'sources/',cfg.myfname],'dataEpoched','recfg')

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

mysurflo = load('markus/canonical_surf_MNI_413.mat');

figure(1); clf
    plot3(mysurflo.surf.pnt(:,1),mysurflo.surf.pnt(:,2),mysurflo.surf.pnt(:,2),'k.')

%% Baseline correction
% When projecting sensor activity into the brain space, deeper sources get 
% larger activity because they're further away from.the sensor (so the 
% dipole-moment has to be much larger to explain a sensor measured field).
% Hence you want to baseline correct and calculate a ratio, or better a 
% statistical change over trials ;-)

% bscfg = [];
% bscfg.demean            = 'yes';        % 'no' or 'yes', whether to apply baseline correction (default = 'no')
% bscfg.baselinewindow    = [-0.2 -0.1];  % [begin end] in seconds, the default is the complete trial (default = 'all')
% [dataEpoched] = ft_preprocessing(bscfg, dataEpoched);

tbs = [-0.2 -0.1]; bs = nan(length(dataEpoched.label),length(dataEpoched.trialinfo));
for tr=1:length(dataEpoched.trialinfo)
    indbs = ( dataEpoched.time{tr} > tbs(1) & dataEpoched.time{tr} < tbs(2) );
    bs(:,tr) = median(dataEpoched.trial{tr}(:,indbs),2) ;
    dataEpoched.trial{tr} = dataEpoched.trial{tr}./repmat(bs(:,tr),1,length(dataEpoched.time{tr}));
end

%% Averaging
clear timelock;
for i=1:2
    tmlkcfg = [];
    tmlkcfg.trials      = ([recfg.megevts.fixStimType] == i & [recfg.megevts.fixType]~=3);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'all';
    tmlkcfg.lpfilter    = 'no';
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'no';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataEpoched); 
end  

%%
% tpnts   = 0.0:0.2:2.4; tw = 0.1;
tbs     = [-0.2 -0.1]; 
for i=1:2
    bs = nan(length(timelock(i).label),length(timelock(i).trialinfo));
    for tr=1:length(timelock(i).trialinfo)
        indbs = ( timelock(i).time > tbs(1) & timelock(i).time < tbs(2) );
        bs(:,tr) = squeeze(median(timelock(i).trial(tr,:,indbs),3)) ;
    end
    avg = nan(length(timelock(i).label),length(timelock(i).time));
    for ch=1:length(timelock(i).label)
        for ts=1:length(timelock(i).time)
            [~,~,~,tmp] = ttest(bs(ch,:)',squeeze(timelock(i).trial(:,ch,ts)));
            timelock(i).avg(ch,ts) = tmp.tstat; 
        end
    end
end


%%
mysurflo.surf.pos = mysurflo.surf.pnt;
gridd.pos       = mysurflo.surf.pnt;

pcfg = [];
pcfg.material   = 'shiny';
pcfg.lightangle = [6,-16];

timepnts = [0.100 0.170]; tw = 0.02;
for i=1:2
    for tp=1:length(timepnts)
        topo = mean(timelock(i).avg(:,timelock(i).time>(timepnts(tp)-tw) & timelock(1).time<(timepnts(tp)+tw)),2);
        plotbraintopo_glossy(mysurflo.surf,topo',[5],pcfg);
        view ([0 0])             % rotate the object in the view
    %     view ([0 -90])             % rotate the object in the view
    end
end

%%    
mysurflo.surf.pos = mysurflo.surf.pnt;
gridd.pos       = mysurflo.surf.pnt;

pcfg = [];
pcfg.material   = 'shiny';
pcfg.lightangle = [6,-16];


% function plotbraintopo(surf,topo,defaultval,mycfg)
% surf must be a 3-d brain model including fields '.pnt',and 'tri'
% topo must be a COLUMN vector (ie. extending horizontally, one row many columns)
% defaultval is an arbitrarily to choose from scaling factor

timepnts = [0.100 0.170]; tw = 0.02;
for i=1%:length(timepnts)
    topo = mean(timelock(1).avg(:,timelock(1).time>(timepnts(i)-tw) & timelock(1).time<(timepnts(i)+tw)),2);
    plotbraintopo_glossy(mysurflo.surf,topo',[20],pcfg);
    view ([0 0])             % rotate the object in the view
%     view ([0 -90])             % rotate the object in the view
end

%%

spcfg = [];
spcfg.method         = 'surface';
spcfg.funparameter   = 'pow';
spcfg.maskparameter  = cfg.funparameter;
spcfg.funcolorlim    = [0.0 1.2];
spcfg.funcolormap    = 'jet';
spcfg.opacitylim     = [0.0 1.2]; 
spcfg.opacitymap     = 'rampup';  
spcfg.projmethod     = 'nearest'; 
spcfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
spcfg.surfdownsample = 10;  % downsample to speed up processing
ft_sourceplot(spcfg, sourceDiffIntNorm);
view ([90 0])             % rotate the object in the view


