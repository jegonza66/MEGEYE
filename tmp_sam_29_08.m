%% Load individual files
% from super_preprocess_simple.m
clear all
close all
clc

cd('D:/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('D:/toolbox/fieldtrip-20180711')

%cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
%addpath('/home/juank/toolbox/fieldtrip-20170827')


ft_defaults;
addpath('my_functions/')
experimentdir = 'D:/Dropbox/big_files_shared/DATA';
directories = {'/13229-001';
               '/13229-001';
               '/13229-001';
               };
sessionfilenames = {'13229-001_MarkusBauer_20180319_04.ds';
                    '13229-001_MarkusBauer_20180319_05.ds';
                    '13229-001_MarkusBauer_20180319_06.ds';
                    };
                
%%
% j = 1
% cfg = [];
%     cfg.directory   = [experimentdir,directories{j}];
% 
%     display('________________________________________________');
%     display(['subject: ',num2str(j)]);
%     cfg.myfname = sessionfilenames{j}([1:(end-3)]);
%     cfg.outname      = [cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat'];
%     
%     cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
%     cfg.continuous      = 'yes';
% 
%     cfg.trialdef.eventtype          = {'UPPT002'};  % Event Channel
%     cfg.trialdef.eventvalue         = [2];          % Values on the event channel
%                                                     % name of the triggers in brain vision format
%     cfg.trialdef.prestim            = 0;   
%     cfg.trialdef.poststim           = 6;
%     cfg = ft_definetrial(cfg);
% 
%     cfg.channel         = {'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
%                                                     % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
%     cfg.bsfilter        = 'no';
%     cfg.demean          = 'no';
%     cfg.detrend         = 'no';
%     data                = ft_preprocessing(cfg);
%     
%     figure(1); clf; plot(data.time{1},data.trial{1}(strcmp(data.label,'UPPT002'),:))

%% First steps of analysis
for j = 1:length(sessionfilenames)
    cfg = [];
    cfg.directory   = [experimentdir,directories{j}];

    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname = sessionfilenames{j}([1:(end-3)]);
    cfg.outname      = [cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat'];
    
    cfg.dataset         = [cfg.directory,filesep,sessionfilenames{j}];
    cfg.continuous      = 'yes';

    cfg.trialdef.eventtype          = {'UPPT002'};  % Event Channel
    cfg.trialdef.eventvalue         = [2];          % Values on the event channel
                                                    % name of the triggers in brain vision format
    cfg.trialdef.prestim            = 0;   
    cfg.trialdef.poststim           = 6;
    cfg = ft_definetrial(cfg);

    cfg.channel         = {'meg';'megref';'UADC*'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
                                                    % (JK, 28/08/2018) The reference is important for proper work of 3rd gradient correction
    cfg.bsfilter        = 'yes';
    cfg.bsfreq          = [49.5,50.5];
    cfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
                                    % (JK, 28/08/2018) Warning: no padding applied because the padding duration is shorter than the trial
    cfg.demean          = 'yes';    % (JK, 22/03/2018) cfg.blc is going to be deprecated
    cfg.detrend         = 'no';
    data                = ft_preprocessing(cfg);
    
    ccfg.resamplefs     = 600;      % (JK, 21/03/2018) I leave it like this, but I guess we'll want to match the eye tracking sampling rate
    ccfg.demean         = 'yes';
    ccfg.detrend        = 'no';     % (JK, 21/03/2018) I leave it like this, but I guess we'll want to change this.
                                    % (JK, 28/08/2018) I change this now.
    data                = ft_resampledata(ccfg, data);

%     % (Markus note) Now the data with reference sensors (aligned with the MEG sensors in time)
%     % and ocular artefact removed are transformed to, e.g. 1st grad
% ??????

    ecfg.channel         = {'UADC*';'UPPT002'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
    edata                = ft_selectdata(ecfg,data);

    mcfg.channel         = {'meg';'megref'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
    data                = ft_selectdata(mcfg,data);
    
    gcfg.gradient       = 'G3BR'; %'none'
    data                = ft_denoise_synthetic(gcfg,data);

    ncfg = [];
    data = ft_appenddata(ncfg, data, edata);    
    
    load('labels_CTF269Nott','labels');
    scfg.channel        = [labels; 'UADC013'; 'UADC014'];
    data                = ft_selectdata(scfg,data);
    save(cfg.outname,'data','-v7.3')
end

%% Merge data files
for j = 1:length(sessionfilenames)
    cfg = [];
    cfg.directory = [experimentdir,directories{j}];

    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname     = sessionfilenames{j}([1:(end-3)]);

    load([cfg.directory,'/matfiles/tmp_eyemap_data_weyech_', cfg.myfname,'.mat']);
    data2merge{j} = data;
end
data = ft_appenddata(cfg, data2merge{1}, data2merge{2}, data2merge{3});
data.fsample = 600; % Add sampling frequency 
clear data2merge

save([cfg.directory,'/matfiles/matias_tmp_eyemap_data_weyech_', cfg.myfname(1:end-3),'_merged.mat'],'data','-v7.3')

%% Change trial markers.
j = 1;
cfg = [];
cfg.directory   = [experimentdir,directories{j}];

display('________________________________________________');
display(['subject: ',num2str(j)]);
cfg.myfname = sessionfilenames{j}([1:(end-6)]);
cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/matias_tmp_eyemap_data_weyech_', cfg.myfname,'_merged.mat'];
cfg.eyemap_sequence         = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];
fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
%% MATIAS 

%cfg = [];
cfg.layout = 'CTF275.lay';
cfg.layout = ft_prepare_layout(cfg);
reading layout from file CTF151.lay
the call to "ft_prepare_layout" took 0 seconds and required the additional allocation of an estimated 1 MB
figure; ft_plot_lay(cfg.layout);


%% Export to EEGLAB
clear data 
addpath('/home/juank/toolbox/eeglab14_1_1b/')
eeglab
EEG = fieldtrip2eeglab(cfg.eyemap_datafname_weyech);
eeglab redraw

%%
load(cfg.eyemap_sequence)
s.info = sequence.info(sequence.value~=200);
s.value = sequence.value(sequence.value~=200);

for i=1:length(s.value)
    EEG.epoch(i).event = i;
    EEG.epoch(i).eventtype = s.value(i);
    EEG.epoch(i).eventlatency = 0;
    EEG.epoch(i).eventurevent = i;
    EEG.epochdescription = sequence.equivalence;

    EEG.event(i).type = s.value(i);
    EEG.event(i).latency = 999999; % This is probably not important, but I want to see an ERROR if it uses it
    EEG.event(i).urevent = i;
    EEG.event(i).epoch = i;
end

EEG = eeg_checkset(EEG);
pop_saveset(EEG,cfg.eyemap_datafname_weyech);
eeglab redraw
% EEG.epoch: 1x30 struct array with fields:
%     event
%     eventtype
%     eventlatency
%     eventurevent

% EEG.event: 1x30 struct array with fields:
%     type
%     latency
%     urevent
%     epoch

%% 
figure(10); clf 
    subplot(1,2,1)
        imagesc(EEG.times,1:30,squeeze(EEG.data(270,:,:))')
    subplot(1,2,2)
        imagesc(EEG.times,1:30,squeeze(EEG.data(271,:,:))')
        
%% http://cognitrn.psych.indiana.edu/busey/temp/eeglabtutorial4.301/maintut/ICA_decomposition.html
% Very important note: We usually run ICA using many more trials that the 
% sample decomposition presented here. As a general rule, finding N stable 
% components (from N-channel data) may require more than 3*N^2 (according to 
% Tony Bell, inventor of ICA Infomax) data sample points (at each channel). 
% In this example with 32 channels, we have 30800 data points and we needed 
% more than 3*32^2=3072 points. Most persons that failed to find reliable 
% component using ICA did not follow this rule.

fprintf('Number of requiered data sample-points: %d\n',3*EEG.nbchan^2)
fprintf('Number of data sample-points that we have: %d\n',length(EEG.times)*EEG.trials)
    