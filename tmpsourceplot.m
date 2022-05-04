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
sessions            = {[4:6]};
isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

%%
isess = 1;
headmodelfile = [cfg.directory,'tmpcovariance',filesep,'fwdmdlspm12raw_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
load(headmodelfile);
positions = forward.forward.mesh.vert;

timelockfile = [cfg.directory,'tmpcovariance',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
load(timelockfile);

leadfieldfile = [cfg.directory,'tmpcovariance',filesep,'leadfieldBNDEM_',num2str(sessions{isubject}(isess)),'.mat'];
load(leadfieldfile); % It has 'grid'

filterfile = [cfg.directory,'tmpcovariance',filesep,'LCMVBNDEMnofiltSAMhires_0',num2str(sessions{isubject}(isess)),'.mat'];


% UNFILTERED LCMV BEAMFORMER
timelock.grad = ft_convert_units(timelock.grad,'m');
srcfg = [];
srcfg.method        = 'sam';
srcfg.feedback      = 'none';
srcfg.channel       = {'MEG'}; %'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
srcfg.vol           = forward.forward.vol;
srcfg.meansphereorigin = nanmean(srcfg.vol.bnd.pos,1);
srcfg.grad          = timelock.grad;
srcfg.parallel      = 'no';
srcfg.keepfilter    = 'yes';
srcfg.fixedori      = 'spinning';
srcfg.grid          = grid;
srcfg.keeptrials    = 'yes';

figure(1); clf
    plot3(srcfg.grid.pos(:,1),srcfg.grid.pos(:,2),srcfg.grid.pos(:,3),'b.')
    hold on; plot3(srcfg.vol.bnd.pos(:,1),srcfg.vol.bnd.pos(:,2),srcfg.vol.bnd.pos(:,3),'r.')
    hold on; plot3(srcfg.grad.chanpos(:,1),srcfg.grad.chanpos(:,2),srcfg.grad.chanpos(:,3),'ko')

if (length(size(timelock.cov))==2)
    srcfg.lambda        = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
elseif (length(size(timelock.cov))==3)
    srcfg.lambda        = nanmean(abs(diag(squeeze(mean(timelock.cov)))))/50; %nanmean(abs(timelock.cov(:)));
end
source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));

% load(filterfile,'source');













cfg.myfname = ['dataEpoched_',sessionfilenames{isubject},'_src413nofilt_merged.mat'];
load([cfg.directory,'sources/',cfg.myfname],'dataEpoched','recfg')
bscfg = [];
bscfg.demean            = 'yes';        % 'no' or 'yes', whether to apply baseline correction (default = 'no')
bscfg.baselinewindow    = [-0.2 -0.1];  % [begin end] in seconds, the default is the complete trial (default = 'all')
[dataEpoched] = ft_preprocessing(bscfg, dataEpoched);

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

mysurflo = load('markus/canonical_surf_MNI_413.mat');

figure(1); clf
    plot3(mysurflo.surf.pnt(:,1),mysurflo.surf.pnt(:,2),mysurflo.surf.pnt(:,2),'k.')

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

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolorlim    = [0.0 1.2];
cfg.funcolormap    = 'jet';
cfg.opacitylim     = [0.0 1.2]; 
cfg.opacitymap     = 'rampup';  
cfg.projmethod     = 'nearest'; 
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfdownsample = 10;  % downsample to speed up processing
ft_sourceplot(cfg, source);
view ([90 0])             % rotate the object in the view