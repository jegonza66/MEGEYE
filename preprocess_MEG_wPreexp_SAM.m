%% Load individual files
% from super_preprocess_simple.m
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath(genpath('/home/juank/toolbox/fieldtrip-20170827'),'-end')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
% experimentdir = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
experimentdir = '/home/juank/Dropbox/big_files_shared/DATA';
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
        fprintf('%s\ņ',cfg.myfname)
        cfg.eyemap_datafname        = [cfg.directory,'/matfiles/tmp_eyemap_data_', cfg.myfname,'.mat'];
        cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat'];
        cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/tmp_eyemap_cfg_',  cfg.myfname,'.mat'];        
        
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
filetypes = {'tmp_eyemap_data_','tmp_eyemap_data_weyech_'};
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
    save([cfg.directory,'/matfiles/',filetypes{i}, cfg.myfname(1:end-3),'_merged.mat'],'data','-v7.3')
end

%% Merge cfg files
j = 1;
directory               = [experimentdir,directories{j}];
myfname                 = sessionfilenames{j}([1:(end-3)]);
eyemap_preprocfname     = [directory,'/matfiles/tmp_eyemap_cfg_',myfname,'.mat'];

myfname_merged          = [sessionfilenames{j}([1:(end-6)]) '_merged'];
dataset                 = [];

load(preprocfname);
cfg.dataset     = dataset;
cfg.datafile    = [directory,'/matfiles/tmp_eyemap_data_',myfname_merged,'.mat'];
cfg.headerfile  = [directory,'/matfiles/tmp_eyemap_data_',myfname_merged,'.mat'];
save([directory,'/matfiles/tmp_eyemap_cfg_',myfname_merged,'.mat'],'cfg');

%% Load merged files
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
% experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
experimentdir = '/home/juank/Dropbox/big_files_shared/DATA';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319_merged.ds'};
                
%% Change trial markers.
j = 1;
cfg.directory = [experimentdir,directories{j}];
cfg.myfname = sessionfilenames{j}([1:(end-3)]);

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);
cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
cfg.eyemap_datafname        = [cfg.directory,'/matfiles/tmp_eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/tmp_eyemap_cfg_',  cfg.myfname,'.mat'];
    
fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
    
%% Export to EEGLAB
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
% experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
experimentdir = '/home/juank/Dropbox/big_files_shared/DATA';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319_merged.ds'};

j = 1;
cfg.directory = [experimentdir,directories{j}];

cfg.myfname = [sessionfilenames{j}(1:(end-6)) '_merged'];

display('________________________________________________');
display(['subject: ',num2str(j) ', ' cfg.myfname]);

cfg.eyemap_datafname        = [cfg.directory,'/matfiles/tmp_eyemap_data_', cfg.myfname,'.mat'];
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat'];
cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/tmp_eyemap_cfg_',  cfg.myfname,'.mat'];

load(cfg.eyemap_datafname);
data.fsample = 600;
save(cfg.eyemap_datafname,'data','-v7.3');
clear data

%%
addpath('/home/juank/toolbox/eeglab14_1_1b/')
eeglab
EEG = fieldtrip2eeglab(cfg.eyemap_datafname);

