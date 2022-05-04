%% Load individual files
% from super_preprocess_simple.m
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath(genpath('/home/juank/toolbox/fieldtrip-20170827'),'-end')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
% addpath('my_functions/')
% experimentdir = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
experimentdir = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories = {'/13229-001';
               '/13229-001';
               '/13229-001';
               };
sessionfilenames = {'13229-001_MarkusBauer_20180319_04.ds';
                    '13229-001_MarkusBauer_20180319_05.ds';
                    '13229-001_MarkusBauer_20180319_06.ds';
                    };
                
%% First steps of analysis
for j = 1:length(sessionfilenames)
% j = 1;
    cfg.directory = [experimentdir,directories{j}];
    
    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname = sessionfilenames{j}([1:(end-3)]);
    
    cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
    
    if exist(cfg.dataset,'dir')
        cfg.datafname       = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
        cfg.datafname_weyech= [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
        cfg.preprocfname    = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];
        
        cfg.eyemap_datafname       = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
        cfg.eyemap_datafname_weyech= [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
        cfg.eyemap_preprocfname    = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];
        
        if ~exist([cfg.directory,'/matfiles'],'dir')
            mkdir([cfg.directory,'/matfiles']);
        end;
    
        % DO THE EPOCHING & BLINK-DETECTION
        % Setup eventmarks (saves a configuration file).
        fun_epoching(cfg);
        
        % Read data and perform some preprocessing (saves a data file).
        fun_readdata_simple(cfg);
    end;
end;

%% Merge data files
filetypes = {'data_','data_weyech_','eyemap_data_','eyemap_data_weyech_'};
for i = 1:length(filetypes)
    for j = 1:length(sessionfilenames)
        cfg.directory = [experimentdir,directories{j}];

        display('________________________________________________');
        display(['subject: ',num2str(j)]);
        cfg.myfname     = sessionfilenames{j}([1:(end-3)]);
        
        load([cfg.directory,'/matfiles/',filetypes{i}, cfg.myfname,'.mat']);
        data2merge{j} = data;
    end
    data = ft_appenddata(cfg, data2merge{1}, data2merge{2}, data2merge{3});
    save([cfg.directory,'/matfiles/',filetypes{i}, cfg.myfname(1:end-3),'_merged.mat'],'data')
end

%% Merge cfg files
j = 1;
directory               = [experimentdir,directories{j}];
myfname                 = sessionfilenames{j}([1:(end-3)]);
preprocfname            = [directory,'/matfiles/cfg_',myfname,'.mat'];
eyemap_preprocfname     = [directory,'/matfiles/eyemap_cfg_',myfname,'.mat'];

myfname_merged          = [sessionfilenames{j}([1:(end-6)]) '_merged'];
dataset                 = [];


load(preprocfname);
cfg.dataset     = dataset;
cfg.datafile    = [directory,'/matfiles/data_',myfname_merged,'.mat'];
cfg.headerfile  = [directory,'/matfiles/data_',myfname_merged,'.mat'];
save([directory,'/matfiles/cfg_', myfname_merged,'.mat'],'cfg');

load(preprocfname);
cfg.dataset     = dataset;
cfg.datafile    = [directory,'/matfiles/eyemap_data_',myfname_merged,'.mat'];
cfg.headerfile  = [directory,'/matfiles/eyemap_data_',myfname_merged,'.mat'];
save([directory,'/matfiles/eyemap_cfg_',myfname_merged,'.mat'],'cfg');

%% Load merged files
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319_merged.ds'};
                
%% Change trial markers.
j = 1;
cfg.directory = [experimentdir,directories{j}];
cfg.myfname = sessionfilenames{j}([1:(end-3)]);

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.sequence                = [cfg.directory,'/matfiles/sequence_', cfg.myfname,'.mat'];
cfg.eyemap_sequence         = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];
    
fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.datafname,                 cfg.datafname,              cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.datafname_weyech,          cfg.datafname_weyech,       cfg.sequence)
    
%% Simple Artifact removal
% from super_preprocess_simple.m
                
j = 1;
cfg.directory = [experimentdir,directories{j}];

cfg.myfname = sessionfilenames{j}([1:(end-3)]);

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.dataset                 = [cfg.directory,filesep,sessionfilenames{j}];
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.datafname_clean         = [cfg.directory,'/matfiles/data_clean_', cfg.myfname,'.mat'];
% cfg.datafname_weyech_clean  = [cfg.directory,'/matfiles/data_weyech_clean_', cfg.myfname,'.mat'];

cfg.eyemap_datafname_clean  = [cfg.directory,'/matfiles/eyemap_data_clean_', cfg.myfname,'.mat'];
% cfg.eyemap_datafname_weyech_clean= [cfg.directory,'/matfiles/eyemap_data_weyech_clean_', cfg.myfname,'.mat'];

% cfgfieldsin     = {'datafname','datafname_weyech','eyemap_datafname','eyemap_datafname_weyech'};
% cfgfieldsout    = {'datafname_clean','datafname_weyech_clean','eyemap_datafname_clean','eyemap_datafname_weyech_clean'};
cfgfieldsin     = {'datafname','eyemap_datafname'};
cfgfieldsout    = {'datafname_clean','eyemap_datafname_clean'};
for i = 1:length(cfgfieldsin)
    % ALLOW REJECTION OF BAD COMPONENTS OR BAD CHANNELS WITH PCA (WHOLE TEMPORARY DATA SET)
    fun_rawdatacomponentrejection(cfg.(cfgfieldsin{i}),cfg.(cfgfieldsout{i}));
end

% cfg.datafname_clean_wEM = 'datafname_clean_wEM';
% fun_rawdatacomponentrejection_usingEyeMap(cfg.datafname,cfg.eyemap_datafname_clean,cfg.datafname_clean_wEM);

%% Less Simple Artifact removal
% I'm going to estimated the PCs in the EyeMap data, and use them to clean
% the Task data                
j = 1;
cfg.directory = [experimentdir,directories{j}];

cfg.myfname = sessionfilenames{j}([1:(end-3)]);

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.dataset                 = [cfg.directory,filesep,sessionfilenames{j}];
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.datafname_clean         = [cfg.directory,'/matfiles/data_clean_', cfg.myfname,'.mat'];
% cfg.datafname_weyech_clean  = [cfg.directory,'/matfiles/data_weyech_clean_', cfg.myfname,'.mat'];

cfg.eyemap_datafname_clean  = [cfg.directory,'/matfiles/eyemap_data_clean_', cfg.myfname,'.mat'];
% cfg.eyemap_datafname_weyech_clean= [cfg.directory,'/matfiles/eyemap_data_weyech_clean_', cfg.myfname,'.mat'];

% ALLOW REJECTION OF BAD COMPONENTS OR BAD CHANNELS WITH PCA (WHOLE TEMPORARY DATA SET)
% [comp,compno]=fun_rawdatacomponentrejection_v3(cfg.eyemap_datafname,cfg.eyemap_datafname_clean,'pca');
[comp,compno]=fun_rawdatacomponentrejection_v3(cfg.eyemap_datafname,cfg.eyemap_datafname_clean,'ica',{100});

load(cfg.eyemap_datafname_weyech)
eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/'; 
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata_EyeMap.mat']);
eyedata = eyedata_EyeMap;

eyedata = fun_estimating_lag(data,eyedata,0);
fprintf('Total fixations = %d\n',sum([eyedata.Nfix]))

recfg = [];
recfg.filterNonStimFix = 0;
recfg.epochLimits_ms= [-500 500]; % In ms,
recfg.thrDuration   = [50 NaN]; % In ms, [200 NaN] means only larger than 200ms
% recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[compEpoched,recfg] = fun_reepoch_data(recfg,comp,eyedata);

%% Less Simple Artifact removal: Appling Plochl et al (2012) criteria
t = compEpoched.time{1};
prevSaccDuration    = [recfg.megevts.prevSaccDuration]/1000;
fixDuration         = [recfg.megevts.fixDuration]/1000;
var_fix     = nan(length(compEpoched.label),length(compEpoched.trial));
var_sacc    = nan(length(compEpoched.label),length(compEpoched.trial));
for ic=1:length(compEpoched.label)
    for tr=1:length(compEpoched.trial)
       var_fix(ic,tr)   = var(compEpoched.trial{tr}( ic , t>0 & t<min([fixDuration(tr),0.5]) ) );
       var_sacc(ic,tr)  = var(compEpoched.trial{tr}( ic , t>(-prevSaccDuration(tr)-0.005) & t<0.010 ) );
    end
end
ratio   = mean(var_sacc./var_fix,2);

thr     = 0.1:0.05:2;
cratio  = nan(size(thr));
for i=1:length(thr); 
    cratio(i) = sum(ratio>thr(i))/length(ratio);
end
figure(1); clf; set(gcf,'Color','w') 
    hold on
        plot([1.1 1.1],[0 1],'k--')
        plot(thr,cratio,'k-','LineWidth',3)
    hold off
    box on
    xlabel('Threshold')
    ylabel('P(var_{sacc}/var_{fix} > Threshold)')

%% Less Simple Artifact removal: Plot partial results
pltcfg.layout   = 'CTF269Nott.lay';
pltcfg.marker   = 'off';
pltcfg.comment  = 'no';
pltcfg.title    = 'off';

[~,ind] = sortrows(ratio,-1);
theta       = [recfg.megevts.prevSaccDirection];
horiz_r     = ( (theta>(2*pi-pi/6)) | (theta<pi/6) );
horiz_l     = ( theta>(pi-pi/6) & theta<(pi+pi/6) );
vert_up     = ( theta>(pi/2 - pi/6) & theta<(pi/2 + pi/6) );
vert_do     = ( theta>(3*pi/2-pi/6) & theta<(3*pi/2+pi/6) );
YLIMI = [-2e-12 2e-12];
icmax = 8;
figure(2); clf; set(gcf,'Color','w','Position',[66 1 1855 1001])
    for ic = 1:icmax
        subplot(icmax,5,5*(ic-1) + 1)
            pltcfg.component = ind(ic);
            ft_topoplotIC(pltcfg, comp)
        
        subplot(icmax,5,5*(ic-1) + 2)
            m=fun_meancurve_cells(compEpoched.trial(vert_up),ind(ic));
            mSaccDur = mean(prevSaccDuration(vert_up));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Upward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); ylabel('Amplitude (fT)'); end
        subplot(icmax,5,5*(ic-1) + 3)
            m=fun_meancurve_cells(compEpoched.trial(vert_do),ind(ic));
            mSaccDur = mean(prevSaccDuration(vert_do));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Downward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
        subplot(icmax,5,5*(ic-1) + 4)
            m=fun_meancurve_cells(compEpoched.trial(horiz_l),ind(ic));
            mSaccDur = mean(prevSaccDuration(horiz_l));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Leftward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
        subplot(icmax,5,5*(ic-1) + 5)
            m=fun_meancurve_cells(compEpoched.trial(horiz_r),ind(ic));
            mSaccDur = mean(prevSaccDuration(horiz_r));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Rightward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
    end
    
%% Less Simple Artifact removal: Checking the blinks    
blicfg = [];
blicfg.filterNonStim = 0;
blicfg.epochLimits_ms= [-1000 1000]; % In ms,
blicfg.thrDuration   = [50 NaN]; % In ms, [200 NaN] means only larger than 200ms
% recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[compEpochedbli,blicfg] = fun_reepoch_data_blinks(blicfg,comp,eyedata);

pltcfg.layout   = 'CTF269Nott.lay';
pltcfg.marker   = 'off';
pltcfg.comment  = 'no';
pltcfg.title    = 'off';

YLIMI = [-2e-12 2e-12];
icmax = 8;
figure(3); clf; set(gcf,'Color','w','Position',[66 1 1855 1001])
    for ic = 1:icmax*3
        subplot(icmax,6,2*(ic-1) + 1)
            pltcfg.component = ic;
            ft_topoplotIC(pltcfg, comp)
        
        subplot(icmax,6,2*(ic-1) + 2)
            m=fun_meancurve_cells(compEpoched.trial,ic);
            hold on
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Blink'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); ylabel('Amplitude (fT)'); end
    end

%% Less Simple Artifact removal: Removing components
compno = find(ratio>1.1);
fun_rawdatacomponentrejection_v3(cfg.eyemap_datafname,cfg.eyemap_datafname_clean,'previous_unmixing',{comp,compno});
fun_rawdatacomponentrejection_v3(cfg.datafname,cfg.datafname_clean,'previous_unmixing',{comp,compno});

% cfg.datafname_clean_wEM = 'datafname_clean_wEM';
% fun_rawdatacomponentrejection_usingEyeMap(cfg.datafname,cfg.eyemap_datafname_clean,cfg.datafname_clean_wEM);

%% Plot some trial to check the simple artifact removal
j = 1;
cfg.directory = [experimentdir,directories{j}];
cfg.myfname = sessionfilenames{j}([1:(end-3)]);
display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.datafname_clean         = [cfg.directory,'/matfiles/data_clean_', cfg.myfname,'.mat'];
cfg.datafname_weyech_clean  = [cfg.directory,'/matfiles/data_weyech_clean_', cfg.myfname,'.mat'];

cfg.eyemap_datafname_clean  = [cfg.directory,'/matfiles/eyemap_data_clean_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech_clean= [cfg.directory,'/matfiles/eyemap_data_weyech_clean_', cfg.myfname,'.mat'];

old = load(cfg.datafname);
new = load(cfg.datafname_clean);
% old = load(cfg.eyemap_datafname);
% new = load(cfg.eyemap_datafname_clean);

tmlkcfg = [];
tmlkcfg.vartrllength = 0;
tmlkcfg.channel = 'MEG';
tmlkcfg.lpfilter = 'yes';
tmlkcfg.lpfreq = 45;
tmlkcfg.keeptrials = 'yes';
tmlkcfg.covariance = 'no';
otimelock = ft_timelockanalysis(tmlkcfg,old.data); 
ntimelock = ft_timelockanalysis(tmlkcfg,new.data); 
% save(timelockfile,'timelock');

indch = sortrows([find(cellfun(@(x) ~isempty(strfind(x,'MLF1')),otimelock.label));...
                find(cellfun(@(x) ~isempty(strfind(x,'MRF1')),otimelock.label))]);
            
tr = 1;
figure(1);clf
    set(gcf,'Color','w')
    axes('Position',[0.075 0.550 0.400 0.400])
        plot(otimelock.time,squeeze(otimelock.trial(tr,indch,:)))
        set(gca,'XLim',[otimelock.time(1) otimelock.time(end)],'XTickLabel',{})
        set(gca,'YLim',[-2*10^-12 2*10^-12])
        ylabel('Amp')
        
    axes('Position',[0.075 0.100 0.400 0.400])
        plot(ntimelock.time,squeeze(ntimelock.trial(tr,indch,:)))
        set(gca,'XLim',[otimelock.time(1) otimelock.time(end)])
        set(gca,'YLim',[-2*10^-12 2*10^-12])
        ylabel('Amp')
        xlabel('Time (S)')

    axes('Position',[0.525 0.550 0.400 0.400])
        imagesc(otimelock.time,1:length(otimelock.label),squeeze(otimelock.trial(tr,:,:)))
        set(gca,'XLim',[otimelock.time(1) otimelock.time(end)],'XTickLabel',{})
        set(gca,'CLim',[-1.5*10^-12 1.5*10^-12])
        set(gca,'YAxisLocation','right','YTickLabel',{})
        ylabel('Channels')
    axes('Position',[0.525 0.100 0.400 0.400])
        imagesc(ntimelock.time,1:length(ntimelock.label),squeeze(ntimelock.trial(tr,:,:)))
        set(gca,'XLim',[otimelock.time(1) otimelock.time(end)])
        set(gca,'CLim',[-1.5*10^-12 1.5*10^-12])
        set(gca,'YAxisLocation','right','YTickLabel',{})
        ylabel('Channels')
        xlabel('Time (S)')

%% Synchronize the MEG time with the ET time
% I moved it to fun_estimating_lag.m
% Detect fixations (I'm going to use EyeLink saccade detection for now).
j = 1;
cfg.directory   = [experimentdir,directories{j}];
cfg.myfname     = sessionfilenames{j}([1:(end-3)]);
cfg.dataset     = [cfg.directory,filesep,sessionfilenames{j}];
% clear experimentdir directories sessionfilenames

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.datafname_clean         = [cfg.directory,'/matfiles/data_clean_', cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.sequence                = [cfg.directory,'/matfiles/sequence_', cfg.myfname,'.mat'];
cfg.eyemap_sequence         = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];

% Load eyemap data
% load(cfg.eyemap_datafname_weyech)
% eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
% subjid = '13229_001';
% load([eyepath subjid '/' subjid '_eyedata_EyeMap.mat']);
% eyedata = eyedata_EyeMap;

% Load task data
data = fun_include_eyech(cfg.datafname_clean,cfg.datafname_weyech);
% load(cfg.datafname_weyech)

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

eyedata = fun_estimating_lag(data,eyedata,0);

recfg = [];
recfg.filterNonStimFix= 1; % In ms,
recfg.epochLimits_ms= [-200 250]; % In ms,
recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataMEG,recfg] = fun_reepoch_data(recfg,data,eyedata);

%% Some analysis: Select channels, Baseline and Average reference.
load('labels_CTF269Nott','labels');
scfg.channel = labels;
dataMEG = ft_selectdata(scfg,dataMEG);

bscfg = [];
% Baseline correction
bscfg.demean            = 'yes';        % 'no' or 'yes', whether to apply baseline correction (default = 'no')
bscfg.baselinewindow    = [-0.2 -0.1];  % [begin end] in seconds, the default is the complete trial (default = 'all')
% % Average reference
% bscfg.reref               = 'yes';        % 'no' or 'yes' (default = 'no')
% bscfg.refchannel          = 'all';        % cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
% bscfg.refmethod           = 'avg';        % 'avg' or 'median' (default = 'avg')

[dataMEG] = ft_preprocessing(bscfg, dataMEG);

%% Some analysis: Select conditions and build the fERFs.
% Faces vs Objects
clear timelock;
for i=1:2
    tmlkcfg = [];
    tmlkcfg.trials = ([recfg.megevts.fixStimType] == i & [recfg.megevts.fixType]~=3);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'MEG';
    tmlkcfg.lpfilter    = 'yes';
    tmlkcfg.lpfreq      = 45;
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'yes';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataMEG); 
end    

%% Some analysis: Plot fERFs imagesc and topographies
% lycfg = [];
% lycfg.layout = 'CTF269Nott.lay';   
% ft_layoutplot(lycfg, data)
label_cond = {'Faces','Objects'};
figure(10); clf
    set(gcf,'Position',[80 1 800 970],'Color','w')
    for i=1:2
        subplot(3,4,4*(i-1)+[1:2]);
            hold on
                title(label_cond{i})
                imagesc(timelock(i).time,1:length(timelock(i).label),timelock(i).avg,[-1.0 1.0]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            set(gca,'XTickLabel',[])
            ylabel('Sensors')
            
        subplot(3,4,4*(i-1)+3);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
        subplot(3,4,4*(i-1)+4);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
    end
        dtimelock = timelock(1);
        dtimelock.avg = timelock(1).avg-timelock(2).avg;
        subplot(3,4,9:10);
            hold on
                title('Faces-Objects')
                imagesc(timelock(1).time,1:length(timelock(1).label),dtimelock.avg,[-1.0 1.0]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            ylabel('Sensors')
            xlabel('Time (s)')

        subplot(3,4,11);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
        subplot(3,4,12);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
            
%% Some analysis: Plot fERFs waveforms          
label_cond = {'Faces','Objects'};
col_cond    = {[1.0 0.0 0.0],[0.0 0.0 1.0]};
shadow_cond = {[1.0 0.5 0.5],[0.5 0.5 1.0]};

indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLO')),data.label)), ...
            find(cellfun(@(x) ~isempty(strfind(x,'MRO')),data.label))}; 
label_chan = {'Left Occ','Right Occ'};
% indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLF1')),data.label)), ...
%             find(cellfun(@(x) ~isempty(strfind(x,'MRF1')),data.label))}; 
% label_chan = {'Left Frontal','Right Frontal'};
figure(11); clf
    set(gcf,'Position',[80 1 800 400],'Color','w')
    
    YLIMI = [-10^-13 10^-13];
    for j=1:2
        subplot(1,2,j);
            hold on
                title(label_chan{j});
%                 plot(timelock(i).time,mean(timelock(1).avg(indch{j},:)) - mean(timelock(2).avg(indch{j},:)),...
%                     'Color','k','LineWidth',2)
                hp = nan(2,1);
                for i=1:2
                    [m,e]=fun_meancurve(timelock(i),indch{j});
                    errorPatch = niceBars2(timelock(i).time,m,e,shadow_cond{i},0.5);hold on;
                    hp(i)=plot(timelock(i).time,m,'Color',col_cond{i},'LineWidth',3);
                    set(errorPatch,'EdgeAlpha',0)
                    
                end
                plot([0 0],YLIMI,'k--')
                plot([min(timelock(i).time) max(timelock(i).time)],[0 0],'k--')
            hold off
            set(gca,'XLim',[min(timelock(i).time) max(timelock(i).time)])
            set(gca,'YLim',[YLIMI])
            xlabel('Time (s)')
            ylabel('Amplitude (fT)')
            if (j==2); legend(hp,label_cond,'FontSize',10); legend boxoff; end
    end
   
%% Some analysis: Stats
stcfg = [];
stcfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
stcfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to 
                                 % evaluate the effect at the sample level
stcfg.correctm = 'cluster';
stcfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                 % will be used for thresholding
stcfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                 % permutation distribution. 
stcfg.minnbchan = 2;               % minimum number of neighborhood channels that is 
                                 % required for a selected sample to be included 
                                 % in the clustering algorithm (default=0).
% cfg.neighbours = neighbours;   % see below
stcfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
stcfg.clustertail = 0;
stcfg.alpha = 0.025;               % alpha level of the permutation test
stcfg.numrandomization = 500;      % number of draws from the permutation distribution

design = zeros(1,size(timelock(1).trial,1) + size(timelock(2).trial,1));
design(1,1:size(timelock(1).trial,1)) = 1;
design(1,(size(timelock(1).trial,1)+1):(size(timelock(1).trial,1) + size(timelock(2).trial,1)))= 2;

stcfg.design = design;             % design matrix
stcfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)

cfg_neighb        = [];
cfg_neighb.method = 'distance';         
cfg_neighb.layout        = 'CTF269Nott.lay';
cfg_neighb.channel       = 'all';
neighbours        = ft_prepare_neighbours(cfg_neighb, dataMEG);

stcfg.neighbours    = neighbours;  % the neighbours specify for each sensor with 
                                 % which other sensors it can form clusters
% ft_neighbourplot(cfg, data);

stcfg.channel       = {'MEG'};     % cell-array with selected channel labels
stcfg.latency       = recfg.epochLimits_ms/1000; % time interval over which the experimental 
                                 % conditions must be compared (in seconds)

[stat] = ft_timelockstatistics(stcfg, timelock(1), timelock(2));   

%% Some analysis: Plot Stats
figure(100); clf
    set(gcf,'Position',[80 1 800 970],'Color','w')
    for i=1:2
        timelockz(i) = timelock(i);
        timelockz(i).avg = timelock(i).avg.*(stat.negclusterslabelmat==1);% | stat.posclusterslabelmat==1);
        subplot(3,4,4*(i-1)+[1:2]);
            hold on
                imagesc(timelock(i).time,1:length(timelock(i).label),timelock(i).avg,[-1.5 1.5]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            set(gca,'XTickLabel',[])
            ylabel('Sensors')
            
        subplot(3,4,4*(i-1)+3);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.5 1.5]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            ft_topoplotER(tpcfg,timelock(i)); %colorbar
            
        subplot(3,4,4*(i-1)+4);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.5 1.5]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            ft_topoplotER(tpcfg,timelock(i)); %colorbar
    end
        subplot(3,4,9:10);
            hold on
                imagesc(timelockz(1).time,1:length(timelockz(1).label),timelockz(1).avg-timelockz(2).avg,[-0.5 0.5]*10^-13)
                plot([0 0],[1 length(timelockz(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelockz(i).time(1) timelockz(i).time(end)])
            set(gca,'YLim',[0.5 length(timelockz(i).label)+0.5])
            ylabel('Sensors')
            xlabel('Time (s)')
                      
%% Target vs dustractors ({'VS','EX','ME'})
clear timelock;
for i=1:2
    tmlkcfg = [];
    tmlkcfg.trials = ([recfg.megevts.fixType] == i & [recfg.megevts.fixTrialType] == 3);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'MEG';
    tmlkcfg.lpfilter    = 'yes';
    tmlkcfg.lpfreq      = 45;
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'yes';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataMEG); 
end    
% label_cond = {'Faces','Objects'};
% col_cond    = {[1.0 0.0 0.0],[0.0 0.0 1.0]};
% shadow_cond = {[1.0 0.5 0.5],[0.5 0.5 1.0]};
label_cond = {'Irrelev','Relev','Targets'};
col_cond    = {[1.0 0.0 0.0],[0.0 1.0 0.0],[0.0 0.0 1.0]};
shadow_cond = {[1.0 0.5 0.5],[0.5 1.0 0.5],[0.5 0.5 1.0]};

% indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLO')),data.label)), ...
%             find(cellfun(@(x) ~isempty(strfind(x,'MRO')),data.label))}; 
label_chan = {[],'MLF','MRF',[],...
              'MLT','MLC','MRC','MRT', ...
              [],'MLP','MRP',[],...
              [],'MLO','MRO',[]};
            
figure(12); clf
    set(gcf,'Position',[80 1 800 400],'Color','w')
    
    YLIMI = [-10^-13 10^-13];
    for j=1:length(label_chan)
        if ~isempty(label_chan{j})
            subplot(4,4,j);
                hold on
                    title(label_chan{j});
                    hp = nan(length(timelock),1);
                    for i=1:length(timelock)
                        indch = find(cellfun(@(x) ~isempty(strfind(x,label_chan{j})),data.label));
                        [m,e]=fun_meancurve(timelock(i),indch);
                        errorPatch = niceBars2(timelock(i).time,m,e,shadow_cond{i},0.5);hold on;
                        hp(i)=plot(timelock(i).time,m,'Color',col_cond{i},'LineWidth',3);
                        set(errorPatch,'EdgeAlpha',0)

                    end
                    plot([0 0],YLIMI,'k--')
                    plot([min(timelock(i).time) max(timelock(i).time)],[0 0],'k--')
                hold off
                set(gca,'XLim',[min(timelock(i).time) max(timelock(i).time)])
                set(gca,'YLim',[YLIMI])
                xlabel('Time (s)')
                ylabel('Amplitude (fT)')
%                 if (j==2); legend(hp,label_cond,'FontSize',10); legend boxoff; end
        end
    end