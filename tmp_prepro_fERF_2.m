% (2018-03-26, Markus) Here is the battery for the SPM pipeline.
%% Set up path
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
% addpath('/home/juank/toolbox/spm12/')
addpath('my_functions/')
addpath('markus/')
ctfdir          = '/home/juank/Dropbox/big_files_shared/CTF_DATA';
etdir           = '/home/juank/Dropbox/big_files_shared/ET_DATA';
rootmridir      = [];%'d:\Projects\MSc_AttentionPrediction\MR_anatomicals\';

subjname        = {'1343301'};
directories     = {''};
sessionfilenames= {{'1343301_MarkusBauer_20180629_01';'1343301_MarkusBauer_20180629_02'}};
sessions        = {[1:2]};

load('labels_CTF269Nott','labels');

subjectinfo.directories         = directories;
subjectinfo.sessionfilenames    = sessionfilenames;

%% First steps of analysis (PREPARE DATA FOR COVARIANCE)
% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [ctfdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    isess = 2
%     for isess = 1:length(sessions{isubject})
        cfg.myfname = sessionfilenames{isubject}{isess};
        cfg.dataset = [cfg.directory,cfg.myfname,'.ds'];
    
        if exist(cfg.dataset,'dir')
            cfg.datafname       = [cfg.directory,'/tmp/data_',        cfg.myfname,'.mat'];
            cfg.datafname_weyech= [cfg.directory,'/tmp/data_weyech_', cfg.myfname,'.mat'];
            cfg.preprocfname    = [cfg.directory,'/tmp/cfg_',         cfg.myfname,'.mat'];

            cfg.eyemap_datafname       = [cfg.directory,'/tmp/eyemap_data_', cfg.myfname,'.mat'];
            cfg.eyemap_datafname_weyech= [cfg.directory,'/tmp/eyemap_data_weyech_', cfg.myfname,'.mat'];
            cfg.eyemap_preprocfname    = [cfg.directory,'/tmp/eyemap_cfg_',  cfg.myfname,'.mat'];

            if ~exist([cfg.directory,'tmpcovariance'],'dir')
                mkdir([cfg.directory,'tmpcovariance']);
            end;
            
            if ~exist([cfg.directory,'tmp'],'dir')
                mkdir([cfg.directory,'tmp']);
            end;

            % DO THE EPOCHING & BLINK-DETECTION
            % Setup eventmarks (saves a configuration file).
            fun_epoching(cfg);

            % Read data and perform some preprocessing (saves a data file).
%             fun_readdata_simple(cfg);
            fun_readdata_simple(cfg,{});
        end;
%     end;

%% Probando los canales del sujeto 4
icfg = cfg;

load(icfg.preprocfname,'cfg');

Rcfg = [];
cfg.channel         = {}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
cfg.bsfilter        = 'yes';
cfg.bsfreq          = [49.5,50.5];
cfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
cfg.demean          = 'yes';    % (JK, 22/03/2018) cfg.blc is going to be deprecated
cfg.detrend         = 'no';
data = ft_preprocessing(cfg);

unique(data.hdr.chantype)

%% Voy a mirar lo que no son:
close all
clc
ind = find(~ismember(data.hdr.chantype, {'clock','headloc','headloc_gof','meggrad','refgrad','refmag'}));

figure; k=0;
for i=1:length(ind)
    d = nan(length(data.trial),size(data.trial{1},2)); for j=1:length(data.trial); d(j,:) = data.trial{j}(ind(i),:); end
    
    if sum(d(:)~=0)>0
%         if (mod(i,49)==1); figure; k=1;
%         else; 
            k=k+1;
%         end
%         subplot(7,7,k)
        subplot(3,4,k)
            hold on
                title(data.label{ind(i)})
                imagesc(data.time{1},1:length(data.trial),d);
            hold off
            set(gca,'YTickLabel',[])
            set(gca,'XTick',[0 4 6]);
            grid on
%             set(gca,'XTickLabel',[])
            axis tight

        fprintf('%s: any 0 = %d, mean = %0.4f, std = %0.4f\n',data.label{ind(i)},sum(d(:)~=0),mean(d(:)),mean(d(:)))
    end
end


%% Candidatos
% UADC001, UADC002, UADC013: Eye movements
% UPPT001, UPPT002
close all

ind = [find(ismember(data.label, 'UPPT001')) find(ismember(data.label, 'UPPT002'))];
ch = 1;
Ntr = length(data.trial);

figure;
    d           = nan(Ntr,size(data.trial{1},2)); for j=1:Ntr; d(j,:) = data.trial{j}(ind(ch),:); end
    d(d<.95)    = 0;
    dl          = (d > .90);
    ddl         = diff(dl')';
    subplot(4,2,1)
        hold on
            imagesc(data.time{1},1:Ntr,d);
            YLIMI = [1 Ntr];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
    subplot(4,2,3)
        hold on
            plot(data.time{1},d')
            YLIMI = [-1 50];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
        ylim(YLIMI)

    subplot(4,2,2)
        hold on
            imagesc(data.time{1},1:Ntr,dl);
            YLIMI = [1 Ntr];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
    subplot(4,2,4)
        hold on
            plot(data.time{1},dl')
            YLIMI = [-0.1 1.1];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
        ylim(YLIMI)
    
    subplot(4,2,6)
        hold on
            imagesc(data.time{1}(1:end-1),1:Ntr,ddl');
            YLIMI = [1 Ntr];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
    subplot(4,2,8)
        hold on
            plot(data.time{1}(1:end-1),ddl)
            YLIMI = [-1.1 1.1];
            plot([0 0],YLIMI,'k--')
            plot([4 4],YLIMI,'k--')
            plot([6 6],YLIMI,'k--')
        hold off
        axis tight
        ylim(YLIMI)
        
Yini = {};
Yend = {};
Tini = [];
Rini = [];
Nini = [];
for j=1:Ntr%-1
    y = ddl(j,:);
    Yini = [Yini find(y==1)];
    Yend = [Yend find(y==-1)];
    Nini = [Nini sum(y==1)];
    
    i0=find(y==1,1,'first')+1;
        Tini = [Tini data.time{1}(i0)];
        Rini = [Rini mean(d(j,i0:(i0+10)))];
end
Rini(Rini>0.9 & Rini<1.1)   = 1;
Rini(Rini>31 & Rini<33)     = 32;

[Tini' Rini']
[Tini' Rini' Nini']
    

%% Edit sequences to trim individual parts.
% % for isubject = 1:length(directories)
% isubject = 1;
% cfg.directory = [experimentdir,directories{isubject},filesep];
% 
% display('________________________________________________');
% display(['subject: ',num2str(isubject)]);
%     
% for isess = 1:length(sessions{isubject})
%     cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess))];
% 
%     display('________________________________________________');
%     display(['subject: ' cfg.myfname]);
% 
%     cfg.gralsequence            = [cfg.directory,'/matfiles/sequence_', cfg.myfname(1:(end-3)),'_merged.mat'];
%     cfg.graleyemap_sequence     = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname(1:(end-3)),'_merged.mat'];
% 
%     cfg.sequence            = [cfg.directory,'/tmp/sequence_', cfg.myfname,'.mat'];
%     cfg.eyemap_sequence     = [cfg.directory,'/tmp/eyemap_sequence_', cfg.myfname,'.mat'];
%     
%     % I've to first load the sequences and keep only the ones I wanted.
%     load(cfg.gralsequence)
%         if (isess==1)
%             ind = 1:60;
%         elseif (isess==2)
%             ind = 61:120;
%         elseif (isess==3)
%             ind = 121:160;
%         end
%         sequence.value  = sequence.value(ind);
%         sequence.info   = sequence.info(ind);
%         save(cfg.sequence,'sequence')
%         clear sequence
%         
%     load(cfg.graleyemap_sequence)
%         if (isess==1)
%             ind = 1:20;
%         elseif (isess==2)
%             ind = 21:40;
%         elseif (isess==3)
%             ind = 41:60;
%         end
% 
%         sequence.value  = sequence.value(ind);
%         sequence.info   = sequence.info(ind);
%         save(cfg.eyemap_sequence,'sequence')
%         clear sequence
% end

%% Change trial markers (and remove OFF markers from the eyemap file).
% for isubject = 1:length(directories)
isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);
    
for isess = 1:length(sessions{isubject})
    cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess))];

    display('________________________________________________');
    display(['subject: ' cfg.myfname]);
%     cfg.dataset = [cfg.directory,cfg.myfname,'.ds'];
    cfg.datafname               = [cfg.directory,'/tmp/data_', cfg.myfname,'.mat'];
    cfg.datafname_weyech        = [cfg.directory,'/tmp/data_weyech_', cfg.myfname,'.mat'];
    cfg.preprocfname            = [cfg.directory,'/tmp/cfg_',  cfg.myfname,'.mat'];

%     cfg.eyemap_datafname        = [cfg.directory,'/tmpcovariance/eyemap_data_', cfg.myfname,'.mat'];
%     cfg.eyemap_datafname_weyech = [cfg.directory,'/tmpcovariance/eyemap_data_weyech_', cfg.myfname,'.mat'];
%     cfg.eyemap_preprocfname     = [cfg.directory,'/tmpcovariance/eyemap_cfg_',  cfg.myfname,'.mat'];

    cfg.sequence            = [cfg.directory,'/tmp/sequence_', cfg.myfname,'.mat'];
%     cfg.eyemap_sequence     = [cfg.directory,'/tmpcovariance/eyemap_sequence_', cfg.myfname,'.mat'];
    
%     fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
%     fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
    fun_replace_eventmarkers(cfg.datafname,                 cfg.datafname,              cfg.sequence)
    fun_replace_eventmarkers(cfg.datafname_weyech,          cfg.datafname_weyech,       cfg.sequence)
    
end

%% Remove components
load('lastcomp','comp','compno')
% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];
        fun_rawdatacomponentrejection_v3(datafname,datafname,'previous_unmixing',{comp,compno});
% 
%         cfg.myfname = ['eyemap_data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
%         datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];
%         fun_rawdatacomponentrejection_v3(datafname,datafname,'previous_unmixing',{comp,compno});
    end

%% Remove extra sensors and Low pass
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];
        load(datafname,'data')
        load('labels_CTF269Nott','labels');
        scfg.channel = labels;
        data = ft_selectdata(scfg,data);

        bscfg = [];
        bscfg.lpfiler = 'yes';        
        bscfg.lpfreq  = 40;
        [data] = ft_preprocessing(bscfg, data);
        
        save(datafname,'data')
    end

%% Before following the analysis I'm need to join the three parts of the experiment.
% I plan to do it straightforward, as I did it before, but we shall see
% what happen with the sources.

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);

for isess = 1:length(sessions{isubject})
    cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
    datafname = [cfg.directory,'tmp',filesep,cfg.myfname];

    load(datafname,'data');
    data2merge{isess} = data;
end
data        = ft_appenddata([], data2merge{1}, data2merge{2}, data2merge{3});
cfg.myfname = ['data_',sessionfilenames{isubject},'_merged.mat'];
save([cfg.directory,'/tmp/',cfg.myfname],'data')

%% For the synchronization I need the eye channels, 
% so I'm going to merge here 
isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);

filetypes = {'data_weyech_'};
for isess = 1:length(sessions{isubject})
    cfg.myfname = ['data_weyech_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];

    load([cfg.directory,'/tmp/',cfg.myfname]);
    data2merge{isess} = data;
end
data = ft_appenddata([], data2merge{1}, data2merge{2}, data2merge{3});
save([cfg.directory,'/tmp/',cfg.myfname(1:end-7),'_merged.mat'],'data')

%% Synchronize the MEG time with the ET time
% clear all
% close all
% clc
% 
% cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath('/home/juank/toolbox/fieldtrip-20170827')
% ft_defaults;
% addpath('my_functions/')
% experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
% directories         = {'/13229-001'};
% sessionfilenames    = {'13229-001_MarkusBauer_20180319'};
% 
% isubject = 1;
% cfg.directory = [experimentdir,directories{isubject},filesep];

% Load merged files
cfg.myfname = ['data_weyech_',sessionfilenames{isubject},'_merged.mat'];
load([cfg.directory,'tmp/',cfg.myfname],'data')

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

data.fsample = round(mean(1./diff(data.time{1})));
[eyedata,xc_xy] = fun_estimating_lag(data,eyedata,0);
fprintf('Total fixations = %d\n',sum([eyedata.Nfix]))

% Load task data
cfg.myfname = ['data_',sessionfilenames{isubject},'_merged.mat'];
load([cfg.directory,'tmp/',cfg.myfname],'data')
data.fsample = round(mean(1./diff(data.time{1})));

recfg = [];
recfg.filterNonStimFix = 1;
recfg.epochLimits_ms= [-200 250]; % In ms,
recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

cfg.myfname = ['dataEpoched250_',sessionfilenames{isubject},'_merged.mat'];
save([cfg.directory,'tmp/',cfg.myfname],'dataEpoched','recfg') 

recfg = [];
recfg.filterNonStimFix = 1;
recfg.epochLimits_ms= [-200 300]; % In ms,
recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

cfg.myfname = ['dataEpoched300_',sessionfilenames{isubject},'_merged.mat'];
save([cfg.directory,'tmp/',cfg.myfname],'dataEpoched','recfg') 

recfg = [];
recfg.filterNonStimFix = 1;
recfg.epochLimits_ms= [-200 400]; % In ms,
recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

cfg.myfname = ['dataEpoched400_',sessionfilenames{isubject},'_merged.mat'];
save([cfg.directory,'tmp/',cfg.myfname],'dataEpoched','recfg') 
