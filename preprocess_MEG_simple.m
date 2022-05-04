% from super_preprocess_simple.m
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath(genpath('/home/juank/toolbox/fieldtrip-20170827'),'-end')
% addpath('/home/juank/toolbox/fieldtrip-20170827')
% ft_defaults;
% addpath('my_functions/')
experimentdir = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories = {'';
               '';
               '';
               };
sessionfilenames = {'13229-001_MarkusBauer_20180319_04.ds';
                    '13229-001_MarkusBauer_20180319_05.ds';
                    '13229-001_MarkusBauer_20180319_06.ds';
                    };
                
%%
for j = 1:length(sessionfilenames)
% j = 1;
    cfg.directory = [experimentdir,directories{j}];
    
    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname = sessionfilenames{j}([1:(end-3)]);
    
    cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
    
    if exist(cfg.dataset,'dir')
        cfg.datafname       = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
        cfg.datafname_weyes = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
        cfg.preprocfname    = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'_weyes.mat'];
        if ~exist([cfg.directory,'/matfiles'],'dir')
            mkdir([cfg.directory,'/matfiles']);
        end;
    
% hdr = ft_read_header(cfg.dataset);
        
        % DO THE EPOCHING & BLINK-DETECTION
        % Setup eventmarks (saves a configuration file).
        fun_epoching(cfg);
        
        % Read data and perform some preprocessing (saves a data file).
        fun_readdata(cfg);
        
%         % ALLOW REJECTION OF BAD COMPONENTS OR BAD CHANNELS WITH PCA (WHOLE TEMPORARY DATA SET)
%         fun_rawdatacomponentrejection(cfg);
    end;
end;

%%
for j = 1:length(sessionfilenames)
% j = 1;
    cfg.directory = [experimentdir,directories{j}];
    
    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname     = sessionfilenames{j}([1:(end-3)]);
    cfg.dataset     = [cfg.directory,filesep,sessionfilenames{j}];
    cfg.datafname   = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
    
    load(cfg.datafname)
    data2merge{j} = data;
end
data = ft_appenddata(cfg, data2merge{1}, data2merge{2}, data2merge{3});
save([cfg.directory,'/matfiles/data_', cfg.myfname,'_merge.mat'],'data')

%%

timelockfile = [cfg.directory,'/matfiles/data_',cfg.myfname,'tmlklcmvnofilt.mat'];
tmlkcfg = [];
tmlkcfg.vartrllength = 0;
tmlkcfg.channel = 'MEG';
tmlkcfg.lpfilter = 'yes';
tmlkcfg.lpfreq = 45;
tmlkcfg.keeptrials = 'no';
tmlkcfg.covariance = 'yes';
timelock = ft_timelockanalysis(tmlkcfg,data); 
save(timelockfile,'timelock');

timelockfile = [cfg.directory,'/matfiles/data_',cfg.myfname,'tmlklcmvhifrqfilt.mat'];
tmlkcfg = [];
tmlkcfg.vartrllength = 0;
tmlkcfg.channel = 'MEG';
tmlkcfg.hpfilter = 'yes';
tmlkcfg.hpfreq = 40;
tmlkcfg.keeptrials = 'no';
tmlkcfg.covariance = 'yes';
timelock = ft_timelockanalysis(tmlkcfg,data); 
save(timelockfile,'timelock');