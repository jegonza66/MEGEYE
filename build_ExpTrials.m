%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load, filter, accomodate and merge multiple files in one dataset per participant
%% File path
% from super_preprocess_simple.m
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
ctfdir                      = session_path.raw;
etdir                       = '/home/juank/Dropbox/big_files_shared/ET_DATA_raw';%session_path.rawet;
matdir                      = session_path.matfiles;
% directories                 = session_path.directories;
% sessionfilenames            = session_path.sessionfilenames;

% (24/08/2017) I removed two empty subjects, check!!!
ind = cellfun(@(x) length(x{1})==0,session_path.behavfilenames);
    session_path.subjname(ind) = [];
    session_path.subjcode(ind) = [];
    session_path.sessionfilenames(ind) = [];
    session_path.behavfilenames(ind) = [];
    
%% Step 0.1. Build unique ExpTrial files.
% for su = 5:length(session_path.subjname)
%     subjname            = session_path.subjname{su};
%     subjcode            = session_path.subjcode{su};
%     sessionfilenames    = session_path.sessionfilenames{su};
%     behavfilenames      = session_path.behavfilenames{su};
% 
%     fprintf('Starting with subject %s\n',subjname)
%     % For the main task
%     cfg.matdir          = [matdir,filesep,subjname];
%     cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];
% 
%     if length(behavfilenames)==1
%         eye = load(fullfile(etdir,[behavfilenames{1},'.mat']));
%         fprintf('%s: %d\n',subjname, length(eye.ExpTrials))
%     end
% end
    
%% KS1, KS2, KS3 ok

%% KS4
su = 4;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'ExpTrials');
                        
%% KS5 Ojo, dos carpetas
su = 5;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    
    TrialPath = [etdir,filesep,behavfilenames{2},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{2},'_Trial*.mat']); TrialNames = {tmp.name};
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'ExpTrials');

%% KS6
su = 6;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% KS7
su = 7;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');
    
%% KS8
su = 8;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');
    
%% 11804
su = 9;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');
    
%% 13130
su = 10;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% 13416
su = 11;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% 13555 Ojo, dos carpetas
su = 12;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% 13556
su = 13;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% 13557
su = 14;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');

%% 13575
su = 15;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        load([TrialPath,filesep,TrialNames{tr}])
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');
%% 10925
su = 17;
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

    TrialPath = [etdir,filesep,behavfilenames{1},'_Trials/'];
    tmp=dir([TrialPath,filesep,behavfilenames{1},'_Trial*.mat']); TrialNames = {tmp.name};
    %MJI warning: the order of tr doesn't follow trial number
    itr_real=zeros(length(TrialNames),1);
    [B,index_sorted_itr]=sort(itr_real);
    for tr = 1:length(TrialNames);
        tok = regexp(TrialNames{tr},'(\d+)+(\D+)+(\d+)+(\D+)+(\d+)','tokens');
        itr_real(tr)=str2double(tok{1}{5}); %finds 'true' trial number
     end
    ExpTrials = [];
    for tr = 1:length(TrialNames);
        tr_real=itr_real(index_sorted_itr(tr)); %not used...keep for later?
        load([TrialPath,filesep,TrialNames{index_sorted_itr(tr)}]);
        ExpTrials = [ExpTrials SingleTrial];
    end
    tmp = load(fullfile(etdir,[behavfilenames{1},'.mat']));
    subj = tmp.subj;
    cfg = tmp.cfg;
    el = tmp.el;
    NEW_TRIALS_NEEDED=~isequal(tmp.ExpTrials,ExpTrials); %if 0, then this extra step 
    %wasn't needed as the original structure contained the correct trials already.
    save(fullfile(etdir,[behavfilenames{1},'m.mat']),'subj','cfg','el','ExpTrials');