% Script for preanalysis of MEG data, see GoogleDocs...MEG_analysis
%2018_MEG_codes/MEG_pipeline
% JUAN: MEG_pipeline is NOT a shared doc. Should we use MEG_analysis instead?
% JUAN: session_path.behavfilenames->this variable doesn't follow the
% Master Excel file. It should contain the participant code (SPMIC)
% followed by -xyz the number of experiment as recorded by the SPMIC
% (2019-02-01, JK) I started from a copy of preprocess_MEG_wPreexp_6_multiuser.m
% (2019-02-28, MJI) Incorporates issues found in joint Juan-Matias skype talk on 27/02/19 
%
% DEFINE MAIN VARIABLES
%sequence.value (1=VS, 2=EX, 3=MEM)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load, filter, accomodate and merge multiple files in one dataset per participant
%% File path
% from super_preprocess_simple.m
clear all
close all
clc

[runpath,code_path,session_path,whoisrunning]=add_paths_matlab_MEG(1);
% checks user identity and define paths according to local environment

cd(runpath)
fields = fieldnames(code_path);
for ii=1:numel(fields)
    dirtoadd=code_path.(fields{ii});        % this adds a dynamic field name
%     addpath(genpath(fullfile(dirtoadd))); % there's no need to add eeglab
                                            % or fieldtrip subdirectories anymore
    addpath(fullfile(dirtoadd),'-end');
end
% if any FieldTrip analysis is needed, then we need to execute ft_defaults
% to include all remaining subdirectories (note that we didn't use
% addpath(genpath()).
% ft_defaults;  % to include fieldtrip subdirectories that are needed
% eeglab;       % to include eeglab subdirectories that are needed

ctfdir                      = session_path.raw;
etdir                       = session_path.rawet;
matdir                      = session_path.matfiles;
% directories                 = session_path.directories; %MIJK delete?
% sessionfilenames            = session_path.sessionfilenames; %MIJK delete?

%Empty subjects are removed and info copied to session_path.exclude  
[session_path] = remove_empty_subjects(session_path);
    
%% Step 0. Check Participant MEG, EYE MOVEMENT, BEHAVIOR
% Behavior
% from 'check_button_press_MEG.m
calc_behavior = 1;
if calc_behavior
    for su = 1:length(session_path.subjname)
        fun_check_behavior(su,session_path,runpath)
    end
end

%% Step 0.1 Check the triggers channel (UPPT002)
% xh: MATLAB trigger codes
% And they were saved in each experiment in the cfg structure:
% % EEG (8bits)
% % 1,2,4,8,16,32,64,128
% % 11 Event Marks
%     cfg.evts.BgnBlock       = 64; %161;
%     cfg.evts.EndBlock       = 128; %162;
%     cfg.evts.BgnTrial       = 16; %151;
%     cfg.evts.EndTrial       = 32; %152;
%     cfg.evts.OnTargetPre    = 8;  % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OnFixCross     = 8;  % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OnDisplay      = 1;  % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OnRetain       = 8;  % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OnTargetPost   = 8;  % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OnResponse     = 8; % + {Task type, Stim type, Resp type} = 3*2*2
%     cfg.evts.OffResponse    = 8; % + {Task type, Stim type, Resp type} = 3*2*2
% %     cfg.evts.notes = [  '{FACE:0, HOUSE:1} * 100 + ' ...
% %                         '{VS-Present:1, VS-Absent:2, EX:3, ME-Present:4, ME-Absent:5} * 10 + '...
% %                         '{target presentation:1, fixation dot:2, display:3, retain:4, question:5, '...
% %                         'respond:6, end:7} * 1 \n'...
% %                         'For example, FACE-EX-display = 000+30+3 = KEY 133'];
% %     {11-17, 21-27, 31-37, 41-47, 51-57, 61-67, ...
% %       111-117, 121-127, 131-137, 141-147, 151-157, 161-167}
%     
% % 12 Event Marks
%     cfg.evts.BgnEyeMap      = 2; %201;
%     cfg.evts.EndEyeMap      = 4; %202;
%     cfg.evts.OnHL           = 211;
%     cfg.evts.OffHL          = 212;
%     cfg.evts.OnHS           = 221;
%     cfg.evts.OffHS          = 222;
%     cfg.evts.OnVL           = 231;
%     cfg.evts.OffVL          = 232;
%     cfg.evts.OnVS           = 241;
%     cfg.evts.OffVS          = 242;
%     cfg.evts.OnBL           = 251;
%     cfg.evts.OffBL          = 252;





xh = [0,1,2,4,8,16,32,64,128,211,212,221,222,231,232,241,242,251,252,253]; % 0 and 253 are not important, I use them just to collect data outside the range
YH = [];
ID = [];
for su = 17 %1:length(session_path.subjname)
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    for j = 1:length(sessionfilenames)
        ID = [ID;[su j]];        
        try
            cfg = [];
            cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
            if ~exist(cfg.matdir,'dir')
                mkdir(cfg.matdir);
            end;

            display('________________________________________________');
            display(['subject: ',subjname]);
            cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
            cfg.dataset         = [ctfdir,filesep,sessionfilenames{j}];
            cfg.continuous      = 'yes';
            cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                               % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
            cfg.bsfilter        = 'no';
            cfg.demean          = 'no';
            cfg.detrend         = 'no';
            data                = ft_preprocessing(cfg);

            yh = hist(data.trial{1},xh);
            YH = [YH; yh];
            
            if (sum(yh)==0)
                unique(data.trial{1});
            end            
        catch
            YH = [YH; nan(size(xh))];
        end
%             figure(1); clf; plot(data.time{1},data.trial{1})
    end;
end

XH2 = xh(2:end-1);
YH2 = YH(:,2:end-1);
figure; 
    subplot(2,1,1); imagesc(xh,1:size(YH2,1),YH2)
    subplot(2,2,4); plot(xh(XH2>100),YH2(:,XH2>100),'.')

% figure; 
%     plot(data.time{1},xh(data.trial{1})) %MJI, ??? error

figure;
     plot(data.time{1},data.trial{1})
%% 
%% Triggers for EyeMap sequence 
% MJI, 01.03.19
% Check in the run code whether the sequence below is fixed or randomised
% across subjects. DONE. The sequence is randomised so it needs to be read
% directly from the pulses.
% % 12 Event Marks
%     cfg.evts.BgnEyeMap      = 2; %201;
%     cfg.evts.EndEyeMap      = 4; %202;
%     cfg.evts.OnHL           = 211;
%     cfg.evts.OffHL          = 212;
%     cfg.evts.OnHS           = 221;
%     cfg.evts.OffHS          = 222;
%     cfg.evts.OnVL           = 231;
%     cfg.evts.OffVL          = 232;
%     cfg.evts.OnVS           = 241;
%     cfg.evts.OffVS          = 242;
%     cfg.evts.OnBL           = 251;
%     cfg.evts.OffBL          = 252;

xh = [0,1,2,4,8,16,32,64,128,211,212,221,222,231,232,241,242,251,252,253]; % 0 and 253 are not important, I use them just to collect data outside the range
YH = [];
ID = [];

    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};

    for j = 1:length(sessionfilenames)
        ID = [ID;[su j]];        
        try
            cfg = [];
            cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
            if ~exist(cfg.matdir,'dir')
                mkdir(cfg.matdir);
            end;

            display('________________________________________________');
            display(['subject: ',sudatabjname]);
            cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
            cfg.dataset         = [ctfdir,filesep,sessionfilenames{j}];
            cfg.continuous      = 'yes';
            cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                               % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
            cfg.bsfilter        = 'no';
            cfg.demean          = 'no';
            cfg.detrend         = 'no';
            data                = ft_preprocessing(cfg);

            BgnEyeMap = find(data.trial{1}==2)
            [indr,indc,valuesEyeMap] = find(data.trial{1});
            [allEyeMap, iaallEyeMap icallEyeMap] = unique(valuesEyeMap,'stable');
            
            
            yh = hist(data.trial{1},xh);
            YH = [YH; yh];
            
            if (sum(yh)==0)
                unique(data.trial{1});
            end            
        catch
            YH = [YH; nan(size(xh))];
        end
%             figure(1); clf; plot(data.time{1},data.trial{1})
    end;


cfg.matdir          = [matdir,filesep,subjname];
cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];
%cfg.datafname
cfg.eyemap_sequence         = [cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];

display('________________________________________________');
            display(['subject: ',subjname]);
            cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
            cfg.dataset         = [ctfdir,filesep,sessionfilenames{j}];
            
            cfg.continuous      = 'yes';
            cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                               % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
            cfg.bsfilter        = 'no';
            cfg.demean          = 'no';
            cfg.detrend         = 'no';
            data                = ft_preprocessing(cfg);



if ~exist([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'file')
    sequence.info = {'SH', 'OFF', 'SV', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'BL', 'OFF', ...
        'SH', 'OFF', 'LH', 'OFF', 'BL', 'OFF', 'SV', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF', ...
        'SV', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF'};
    sequence.value = [101, 200, 103, 200, 103, 200, 104, 200, 104, 200, 102, 200, 101, 200, 104, 200, 102, 200, 104, 200, ...
        101, 200, 102, 200, 104, 200, 103, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200, ...
        103, 200, 101, 200, 104, 200, 102, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200]';
    sequence.equivalence = {'SH', 101; 'LH', 102; 'SV', 103; 'LV', 104; 'BL', 105; 'OFF', 200};
    save([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'sequence');
else
    load([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat']);
end




%% Step 0.1. Build unique ExpTrial files.
% build_ExpTrials.m looks for individual trials stored in Trials folder and
% rewrites the ExpTrials structure with those new fields.
% see build_ExpTrials.m    
% MJI: modify and execute build_ExpTrials for last two subjects
    
%% Step 0.2. Define sequence of trials (PSEUDO-MANUALLY).
for su = 17 %[1:3 5:length(session_path.subjname)]
    subjname            = session_path.subjname{su};
    subjcode            = session_path.subjcode{su};
    sessionfilenames    = session_path.sessionfilenames{su};
    behavfilenames      = session_path.behavfilenames{su};
    
    fprintf('Starting with subject %s\n',subjname)
    % For the main task
    cfg.matdir          = [matdir,filesep,subjname];
    cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];
    
    if length(behavfilenames)==1
        %             eye = load(fullfile(etdir,[behavfilenames{1},'.mat'])); %MJI
        %             s17 OK
        eye = load(fullfile(etdir,[behavfilenames{1},'.mat']));
        % arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials);
        
        if ~exist([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'file');
            fprintf('\t%s doesnt exist\n',[cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat']);
            if (length(eye.ExpTrials)==160)
                fprintf('\t%s: %d trials on ExpTrials\n',length(eye.ExpTrials));
                sequence.value = arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials)';
                sequence.equivalence = {eye.cfg.tasks{1}, 1; eye.cfg.tasks{2}, 2;eye.cfg.tasks{3}, 3};
                sequence.info = {};
                for i=1:length(sequence.value)
                    sequence.info{i} = sequence.equivalence{sequence.value(i),1};
                end
                if(~exist([cfg.matdir,'/matfiles/'],'dir'))
                    mkdir([cfg.matdir,'/matfiles/'])
                end
                save([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'sequence');
            elseif (length(dir([TrialPath,filesep,behavfilenames '_Trial*.mat'])) == 160);
                fprintf('\t%s: %d trials on Trials directory\n',length(eye.ExpTrials));
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
                save([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'sequence');
            end
        else
            fprintf('\t%s exists\n',[cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat']);
            %                 eye2 = load([etdir,filesep,behavfilenames,'_v2.mat']);
            %                 load([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat']);
        end
    else
        fprintf('\tERROR: %s has more than one behavior file\n',subjname);
    end
end


%% Step 1.1. Select Participant
su = 17

%% Step 1.2. For ONE participant: First steps of analysis
% Creates cfg (Fieltrip structure). 
% Define trigger (eventvalue:1, DISPLAY ONSET pre:-1sec, post:3sec+1sec) for epoching
% Run preprocessing (read data for all MEG channels+ET channels,
subjname            = session_path.subjname{su};
subjcode            = session_path.subjcode{su};
sessionfilenames    = session_path.sessionfilenames{su};
behavfilenames      = session_path.behavfilenames{su};

cfg = [];
for j = 1:length(sessionfilenames)
    %     cfg.directory   = [experimentdir,directories{j}];
    cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
    if ~exist(cfg.matdir,'dir')
        mkdir(cfg.matdir);
    end
    
    display('________________________________________________');
    display(['subject: ',subjname]);
    cfg.myfname         = sessionfilenames{j}([1:(end-3)]);
    cfg.dataset         = [ctfdir,filesep,sessionfilenames{j}];
    
    if exist(cfg.dataset,'dir')
        cfg.datafname       = [cfg.matdir,'/matfiles/data_', cfg.myfname,'.mat'];
        cfg.datafname_weyech= [cfg.matdir,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
        cfg.datafname_eyech_only=[cfg.matdir,'/matfiles/data_eyech_only_', cfg.myfname,'.mat'];
        cfg.preprocfname    = [cfg.matdir,'/matfiles/cfg_',  cfg.myfname,'.mat'];
        
        cfg.eyemap_datafname       = [cfg.matdir,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
        cfg.eyemap_datafname_weyech= [cfg.matdir,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
        cfg.eyemap_datafname_eyech_only=[cfg.matdir,'/matfiles/eyemap_data_eyech_only_', cfg.myfname,'.mat'];
        cfg.eyemap_preprocfname    = [cfg.matdir,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];
        
        if ~exist([cfg.matdir,'/matfiles'],'dir')
            mkdir([cfg.matdir,'/matfiles']);
        end
        
        % DO THE EPOCHING & BLINK-DETECTION
        % Setup eventmarks (saves a configuration file).
        fun_epoching(cfg);
        
        % Read data and perform some preprocessing (saves a data file).
        fun_readdata_simple(cfg);
    end
end
%% Builds EYEMAP sequence
%cfg.eyemap_preprocfname     = [cfg.matdir,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];
eyemap_cfg = load(cfg.eyemap_preprocfname);
if ~exist([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'file')
    sequence = fun_buildeyemap_sequence(eyemap_cfg.cfg);
    save([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'sequence');
else
    load([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat']);
end

%% Step 1.3. For ONE participant: Merge data files (needed because output from previous steps 
% might contain different files for different consecutive sessions)
    filetypes = {'data_','data_weyech_','data_eyech_only_','eyemap_data_','eyemap_data_weyech_',...
        'eyemap_data_eyech_only_','eyemap_sequence_'};
    for i = 1:length(filetypes)
        cfg = [];
        cfg.matdirmerged = [matdir,filesep,subjname];
        str = 'data = ft_appenddata(cfg';    
        for j = 1:length(sessionfilenames)
            cfg.matdir      = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
    %         cfg.directory = [experimentdir,directories{j}];

            display('________________________________________________');
            display([filetypes{i}(1:end-1) ': subject: ' subjname ', session: ' num2str(j)]);
            cfg.myfname         = sessionfilenames{j}([1:(end-3)]);

            load([cfg.matdir,'/matfiles/',filetypes{i}, cfg.myfname,'.mat']);
            data2merge{j} = data;

            str = [str, ', data2merge{' num2str(j) '}'];
        end
        str = [str, ');'];
        eval(str);
        if ~exist([cfg.matdirmerged,'/matfiles/'],'dir')
            mkdir([cfg.matdirmerged,'/matfiles/']);
        end
        if(exist('sequence','var') && strcmp(filetypes{i},'eyemap_sequence_'))
            save([cfg.matdirmerged,'/matfiles/',filetypes{i}, cfg.myfname(1:end-3),'_merged.mat'],'data','sequence');
        else
            save([cfg.matdirmerged,'/matfiles/',filetypes{i}, cfg.myfname(1:end-3),'_merged.mat'],'data')       
        end
    end

%% Step 1.4. For ONE participant: Merge cfg files
    j = 1;
    directory               = [matdir,filesep,subjname,filesep,sessionfilenames{j}(1:end-3)];
    myfname                 = sessionfilenames{j}([1:(end-3)]);
    preprocfname            = [directory,'/matfiles/cfg_',myfname,'.mat'];
    eyemap_preprocfname     = [directory,'/matfiles/eyemap_cfg_',myfname,'.mat'];

    myfname_merged          = [sessionfilenames{j}([1:(end-6)]) '_merged'];
    dataset                 = [];

 
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
ctfdir                      = session_path.raw;
etdir                       = session_path.rawet;
matdir                      = session_path.matfiles;
% experimentdir               = session_path.raw;
% directories                 = session_path.directories;
% sessionfilenames            = session_path.sessionfilenames;

[session_path] = remove_empty_subjects(session_path);

%% Step 2.1. Select Participant
su = 17;
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

if ~exist([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'file');
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
        save([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat'],'sequence');
    end
else    
        if(~exist([etdir,filesep,behavfilenames,'_v2.mat'],'file'))
            eye2=eye;
        else
            eye2 = load([etdir,filesep,behavfilenames,'_v2.mat']); 
        end
        load([cfg.matdir,'/matfiles/sequence_', cfg.myfname,'.mat']);
end
%% For the EyeMap
% It should be: 
% eye.cfg.em.sttags
%em      = [eye2.ExpTrials.em];
em      = [eye.ExpTrials.em]; %MJIJK

cfg.matdir          = [matdir,filesep,subjname];
cfg.myfname         = [sessionfilenames{1}([1:(end-6)]) '_merged'];

%load([matdir,filesep,subjname,'/matfiles/eyemap_data_weyech_',cfg.myfname,'.mat'])
load([matdir,filesep,subjname,'/matfiles/eyemap_data_eyech_only_',cfg.myfname,'.mat'])

cheogv = find(strcmp(data.label,'UADC001'));
cheogh = find(strcmp(data.label,'UADC002'));
Ntr = length(data.trial);
figure(1); clf
    set(gcf,'Color','w')
    for tr = 1:Ntr
        subplot(Ntr,1,tr)
            tini = em(tr).timing(1:2:end) - em(tr).timing0;
            tend = em(tr).timing(2:2:end) - em(tr).timing0;
            hold on
                plot(data.time{tr},data.trial{tr}(cheogv,:),'r-')
                plot(data.time{tr},data.trial{tr}(cheogh,:),'b-')
                YLIMI = ylim;
%MJI. commented the next 3 lines (out of range) 
%                 for i=1:length(tini)
%                     plot([tini(i) tini(i)],YLIMI,'k-')
%                     plot([tend(i) tend(i)],YLIMI,'k--')
%                 end
            hold off
    end
       
%% Step 2.3. Change trial markers.
    j = 1;
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
%MJI error: cfg.eyemap sequence doesn't exist
    fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
    fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
    fun_replace_eventmarkers(cfg.datafname,                 cfg.datafname,              cfg.sequence)
    fun_replace_eventmarkers(cfg.datafname_weyech,          cfg.datafname_weyech,       cfg.sequence)
%MJI: run without cfg.eyemap_data... lines            
%% Remove extra channels (Only keep the analog eye channels)
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
% reepoch data based on individual fixations
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