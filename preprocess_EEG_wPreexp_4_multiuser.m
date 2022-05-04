%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1. Load, filter, and verify data in one dataset per participant
%  VERSION 06/08/2018

%% Step 1.0. Load individual files
clear all
close all
clc

extractEyeData = 1; %perform EM analysis

[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_EEG(1);
%checks user identity and define paths according to local environment
cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii}); %this adds a dynamic field name
    %     addpath(genpath(fullfile(dirtoadd)));
    addpath(fullfile(dirtoadd));
end
% ft_defaults;
eegdir                      = session_path.raw;
etdir                       = session_path.rawet;
matdir                      = session_path.matfiles;
% directories                 = session_path.directories;
% sessionfilenames            = session_path.sessionfilenames;

%% Step 1.1. Select Participant
su = 4;

%% Step 1.2. For ONE participant: First steps of analysis
subjname            = session_path.subjname{su};
%subjcode            = session_path.subjcode{su};
sessionfilenames    = session_path.sessionfilenames{su};
Trialfilenames      = session_path.Trialfilenames{su};

cfg = [];
%     cfg.directory   = [experimentdir,directories{j}];
cfg.matdir      = [matdir,filesep,subjname];
if ~exist(cfg.matdir,'dir')
    mkdir(cfg.matdir);
end;

display('________________________________________________');
display(['subject: ',subjname]);
cfg.myfname         = subjname;
cfg.dataset         = [eegdir,filesep,sessionfilenames];

if exist(cfg.dataset,'dir')
    cfg.datafname       = fullfile(cfg.matdir,'matfiles',['data_' cfg.myfname '.mat']);
    cfg.datafname_weyech= fullfile(cfg.matdir,'matfiles',['data_weyech_' cfg.myfname '.mat']);
    cfg.preprocfname    = fullfile(cfg.matdir,'matfiles',['cfg_' cfg.myfname '.mat']);
    
    cfg.eyemap_datafname       = fullfile(cfg.matdir,'matfiles',['eyemap_data_' cfg.myfname '.mat']);
    cfg.eyemap_datafname_weyech= fullfile(cfg.matdir,'matfiles',['eyemap_data_weyech_' cfg.myfname '.mat']);
    cfg.eyemap_preprocfname    = fullfile(cfg.matdir,'matfiles',['eyemap_cfg_' cfg.myfname '.mat']);
    
    if ~exist(fullfile(cfg.matdir,'matfiles'),'dir')
        mkdir(fullfile(cfg.matdir,'matfiles'));
    end;
    
    %% Extract eye data events
    whichEyeAnalyze     = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right','right'};
    whichEye            = {'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right' 'right','right'};
    
    if extractEyeData
            % Eye data extraction parameters
            % - whicheye: which eye is going to be the output ('right' or 'left')
            PARAMS.whicheye            = whichEye{su};
            % - whicheyeanalize: which eye analize in the case of binocular data
            % ('right', 'left' or 'both')
            PARAMS.whicheyeanalyze     = whichEyeAnalyze{su};
            [all_eyelink] = read_eyetracker_data_ASCII_edit([sessionfilenames '.asc'],etdir,PARAMS);
            filename=fullfile(cfg.matdir,'matfiles',[sessionfilenames '_all_eyelink.mat']);
            save(filename,'all_eyelink','PARAMS')
    end
    %%
    
    % DO THE EPOCHING & BLINK-DETECTION
    % Setup eventmarks (saves a configuration file).
    fun_epoching(cfg);
    
    % Read data and perform some preprocessing (saves a data file).
    fun_readdata_simple(cfg);
end;

%% Step 1.3. For ONE participant: Merge data files
filetypes = {'data_','data_weyech_','eyemap_data_','eyemap_data_weyech_'};
for i = 2:length(filetypes)
    cfg = [];
    cfg.matdirmerged = [matdir,filesep,subjname];
    str = 'data = ft_appenddata(cfg';
    for j = 1:length(sessionfilenames)
        cfg.matdir      = [matdir,filesep,subjname,filesep,subjname];
        %         cfg.directory = [experimentdir,directories{j}];
        
        display('________________________________________________');
        display([filetypes{i}(1:end-1) ': subject: ' subjname ', session: ' num2str(j)]);
        cfg.myfname         = subjname;
        
        load([cfg.matdir,'/matfiles/',filetypes{i}, cfg.myfname,'.mat']);
        data2merge{j} = data;
        
        str = [str, ', data2merge{' num2str(j) '}'];
    end
    str = [str, ');'];
    eval(str);
    if ~exist([cfg.matdirmerged,'/matfiles/'],'dir')
        mkdir([cfg.matdirmerged,'/matfiles/']);
    end;
    
    save([cfg.matdirmerged,'/matfiles/',filetypes{i}, cfg.myfname(1:end-3),'_merged.mat'],'data')
end

%% Step 1.4. For ONE participant: Merge cfg files
j = 1;
directory               = [matdir,filesep,subjname,filesep,subjname];
myfname                 = subjname;
preprocfname            = [directory,'/matfiles/cfg_',myfname,'.mat'];
eyemap_preprocfname     = [directory,'/matfiles/eyemap_cfg_',myfname,'.mat'];

myfname_merged          = [sessionfilenames{j}([1:(end-6)]) '_merged'];
dataset                 = [];


load(preprocfname);
cfg.dataset     = dataset;
cfg.datafile    = [matdir,filesep,subjname,'/matfiles/data_',myfname_merged,'.mat'];
cfg.headerfile  = [matdir,filesep,subjname,'/matfiles/data_',myfname_merged,'.mat'];
save([matdir,filesep,subjname,'/matfiles/cfg_', myfname_merged,'.mat'],'cfg');

load(preprocfname);
cfg.dataset     = dataset;
cfg.datafile    = [matdir,filesep,subjname,'/matfiles/eyemap_data_',myfname_merged,'.mat'];
cfg.headerfile  = [matdir,filesep,subjname,'/matfiles/eyemap_data_',myfname_merged,'.mat'];
save([matdir,filesep,subjname,'/matfiles/eyemap_cfg_',myfname_merged,'.mat'],'cfg');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2.
%% Step 2.0. Load merged files
clear all
close all
clc

[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG(1);
%checks user identity and define paths according to local environment
cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii}); %this adds a dynamic field name
    %     addpath(genpath(fullfile(dirtoadd)));
    addpath(fullfile(dirtoadd));
end
% ft_defaults;
eegdir                      = session_path.raw;
etdir                       = session_path.rawet;
matdir                      = session_path.matfiles;
% experimentdir               = session_path.raw;
% directories                 = session_path.directories;
% sessionfilenames            = session_path.sessionfilenames;

%% Step 2.1. Select Participant
su = 10;
subjname            = session_path.subjname{su};
subjcode            = session_path.subjcode{su};
sessionfilenames    = session_path.sessionfilenames{su};
behavfilenames      = session_path.behavfilenames{su};
if length(behavfilenames)>1;
    fprintf('WARNING: This participant has the eye movement file splitted\n')
else
    behavfilenames = behavfilenames{1};
end

%% Step 2.2. Define sequence of trials (PSEUDO-MANUALLY).
% For the main task
cfg.matdir          = [matdir,filesep,subjname];
cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

eye = load([etdir,filesep,behavfilenames,'.mat']);
% arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials);

if (length(eye.ExpTrials)==160)
    sequence.value = arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials)';
    sequence.equivalence = {eye.cfg.tasks{1}, 1; eye.cfg.tasks{2}, 2;eye.cfg.tasks{3}, 3};
    sequence.info = {};
    for i=1:length(sequence.value)
        sequence.info{i} = sequence.equivalence{sequence.value(i),1};
    end
    save([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'sequence');
elseif (length(dir([TrialPath,filesep,behavfilenames '_Trial*.mat'])) == 160);
    if ~exist([etdir,filesep,behavfilenames,'_v2.mat'],'file')
        TrialPath = [etdir,filesep,behavfilenames '_Trials/'];
        ExpTrials = [];
        for tr = 1:160;
            TrialName = [TrialPath,filesep,behavfilenames '_Trial',num2str(tr),'.mat'];
            load(TrialName)
            ExpTrials = [ExpTrials SingleTrial];
        end
        save([etdir,filesep,behavfilenames,'_v2.mat'],'ExpTrials');
    else
        load([etdir,filesep,behavfilenames,'_v2.mat']);
    end
    sequence.value = arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),ExpTrials)';
    sequence.equivalence = {eye.cfg.tasks{1}, 1; eye.cfg.tasks{2}, 2;eye.cfg.tasks{3}, 3};
    sequence.info = {};
    for i=1:length(sequence.value)
        sequence.info{i} = sequence.equivalence{sequence.value(i),1};
    end
end
    
% For the EyeMap
cfg.matdir          = [matdir,filesep,subjname];
%     cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];
%
%     load([matdir,filesep,subjname,'/matfiles/eyemap_data_weyech_',cfg.myfname,'.mat'])
%
%     cheogv = find(strcmp(data.label,'UADC001'));
%     cheogh = find(strcmp(data.label,'UADC002'));
%     Ntr = length(data.trial);
%     figure(1); clf
%         set(gcf,'Color','w')
%         for tr = 1:Ntr
%             subplot(Ntr,1,tr)
%                 hold on
%                     plot(data.time{tr},data.trial{tr}(cheogv,:),'r-')
%                     plot(data.time{tr},data.trial{tr}(cheogh,:),'b-')
%                 hold off
%
%         end

%% Step 2.2b. Define sequence of trials (PSEUDO-MANUALLY). There's a problem with the ExpTrials structure, I'm going to inspect the individual Trials.


%% Step 2.3. Change trial markers.
subjname            = session_path.subjname{su};
subjcode            = session_path.subjcode{su};
sessionfilenames    = session_path.sessionfilenames{su};
behavfilenames      = session_path.behavfilenames{su};

cfg.matdir          = [matdir,filesep,subjname];
cfg.myfname         = [sessionfilenames{j}([1:(end-6)]) '_merged'];

display('________________________________________________');
display(['subject: ',   subjname ', ' cfg.myfname]);
%     cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];

cfg.datafname               = [cfg.matdir,'/matfiles/data_',        cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.matdir,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.matdir,'/matfiles/cfg_',         cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.matdir,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.matdir,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.matdir,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.sequence                = [cfg.matdir,'/matfiles/sequence_',    cfg.myfname,'.mat'];
cfg.eyemap_sequence         = [cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];

fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.datafname,                 cfg.datafname,              cfg.sequence)
fun_replace_eventmarkers(cfg.datafname_weyech,          cfg.datafname_weyech,       cfg.sequence)

%% Step 2.4.
j = 1;
cfg.directory = [experimentdir,directories{j}];
cfg.myfname = subjname;

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
fun_replace_eventmarkers(cfg.datafname,                 cfg.datafname,              cfg.sequence)
fun_replace_eventmarkers(cfg.datafname_weyech,          cfg.datafname_weyech,       cfg.sequence)

%% Remove extra channels (Only keep the analog eye channels)
% Detect fixations (I'm going to use EyeLink saccade detection for now).
j = 1;
cfg.directory   = [experimentdir,directories{j}];
cfg.myfname     = subjname;
cfg.dataset     = [cfg.directory,filesep,sessionfilenames{j}];
% clear experimentdir directories sessionfilenames

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];

cfg.sequence                = [cfg.directory,'/matfiles/sequence_', cfg.myfname,'.mat'];
cfg.eyemap_sequence         = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];

data = fun_include_eyech(cfg.datafname,cfg.datafname_weyech);
save(cfg.datafname_weyech,'data')

data = fun_include_eyech(cfg.eyemap_datafname,cfg.eyemap_datafname_weyech);
save(cfg.eyemap_datafname_weyech,'data')

%% Synchronize the MEG time with the ET time
% Load eyemap data
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
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

%%
% Load task data
load(cfg.datafname_weyech)
eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

eyedata = fun_estimating_lag(data,eyedata,0);
fprintf('Total fixations = %d\n',sum([eyedata.Nfix]))

recfg = [];
recfg.filterNonStimFix = 0;
recfg.epochLimits_ms= [-500 500]; % In ms,
recfg.thrDuration   = [50 NaN]; % In ms, [200 NaN] means only larger than 200ms
% recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

%% Simple Artifact removal
% from super_preprocess_simple.m

j = 1;
cfg.directory = [experimentdir,directories{j}];

cfg.myfname = subjname;

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

%% Plot some trial to check the simple artifact removal
j = 1;
cfg.directory = [experimentdir,directories{j}];
cfg.myfname = subjname;
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

% old = load(cfg.datafname);
% new = load(cfg.datafname_clean);

old = load(cfg.eyemap_datafname);
new = load(cfg.eyemap_datafname_clean);

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

%% Some analysis: Plot fERFs
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

%%
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
stcfg.numrandomization = 100;      % number of draws from the permutation distribution

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
    timelockz(i).avg = timelock(i).avg.*(stat.negclusterslabelmat==1 | stat.posclusterslabelmat==1);
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