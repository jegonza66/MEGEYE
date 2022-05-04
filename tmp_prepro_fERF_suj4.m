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
    
    isess = 1
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

%             % DO THE EPOCHING & BLINK-DETECTION
%             % Setup eventmarks (saves a configuration file).
%             fun_epoching(cfg);
% 
%             % Read data and perform some preprocessing (saves a data file).
% %             fun_readdata_simple(cfg);
%             fun_readdata_simple(cfg,{});
        end;
%     end;




%% define trials / epochs
%     cfg.trialdef.eventtype          = 'gui';
icfg.dataset                    = cfg.dataset;
icfg.continuous                 = 'yes';

icfg.trialdef.eventtype          = {'UPPT002'};  % Event Channel
icfg.trialdef.eventvalue         = [1];          % Values on the event channel
                                                % name of the triggers in brain vision format
icfg.trialdef.prestim            = 1.2;   
icfg.trialdef.poststim           = 8;           %4
icfg = ft_definetrial(icfg);


icfg.channel         = {};

% icfg.bsfilter        = 'yes';
% icfg.bsfreq          = [49.5,50.5];
% icfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
icfg.demean          = 'no';    % (JK, 22/03/2018) cfg.blc is going to be deprecated
icfg.detrend         = 'no';
data = ft_preprocessing(icfg);
    
%% Probando los canales del sujeto 4


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
