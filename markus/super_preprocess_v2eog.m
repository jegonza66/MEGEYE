
function super_preprocess_simple_v2eog(mysubjects)

experimentdir = 'd:\Projects\MSc_AttentionPrediction\MEG\';
[directories,sessionfilenames,sessions] = collectsubjectinfo;

if nargin<1
    mysubjects = [1:length(directories)];
end;
mysubjects

epoching = 0;
reading = 0;
rejection = 1;

for isubject = mysubjects
    
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        
        if exist(cfg.dataset,'dir')
            
            cfg.preprocfname = [cfg.directory,'matfiles',filesep,'cfg_',num2str(isess),'.mat'];
            cfg.prepreyefname = [cfg.directory,'matfiles',filesep,'cfgEYE_',num2str(isess),'.mat'];
            
            cfg.FTdatafname = [cfg.directory,'matfiles',filesep,'data_',num2str(isess),'.mat'];
            cfg.FTeyedatafname = [cfg.directory,'matfiles',filesep,'dataEYE_',num2str(isess),'.mat'];
            
            if ~exist([cfg.directory,'matfiles'],'dir')
                mkdir([cfg.directory,'matfiles']);
            end;
            
            if epoching, epoching_blinkdetect(cfg); end;
            if reading, readdata(cfg); end;
            if rejection, rawdatacomponentrejection(cfg); end;
            
        end;
        
    end;
    
end;

end



%% DO THE EPOCHING & BLINK-DETECTION

function epoching_blinkdetect(icfg)

cfg.dataset                    = icfg.dataset;
cfg.continuous                 = 'yes';

%% define trials / epochs
cfg.trialdef.eventvalue         = [11:16,21:26,31:36,41:46,51:56,61:66,71:76]+10; %10 times the probab + stim id (face/house) + coherence level
cfg.trialdef.eventtype          = {'UPPT002'};
cfg.trialdef.prestim            = 1.3;
cfg.trialdef.poststim           = 2;
cfg = ft_definetrial(cfg);
%cfg.trl([1,end],:) = [];
save(icfg.preprocfname,'cfg');

cfg.bsfilter                    = 'yes';
cfg.bsfreq                      = [49,51];
cfg.padding                     = 5;
cfg.blc                         = 'yes';%'no';
cfg.detrend                     = 'no';
cfg.artfctdef.feedback          = 'yes';
cfg.artfctdef.reject            = 'partial';
cfg.artfctdef.minaccepttim      = 0.1;
cfg.artfctdef.type              = {'eog'};
cfg.artfctdef.eog.bpfilter     = 'yes';
cfg.artfctdef.eog.bpfilttype   = 'but';
cfg.artfctdef.eog.bpfreq       = [1 15];
cfg.artfctdef.eog.bpfiltord    = 4;
cfg.artfctdef.eog.hilbert      = 'yes';
cfg.artfctdef.eog.channel      = {'EEG057';'EEG058'};
cfg.artfctdef.eog.cutoff       = 4;   % z-value at which to threshold
cfg.artfctdef.eog.trlpadding   = 0;   % indicates range to search 4 artifact
cfg.artfctdef.eog.fltpadding   = 1;   % indicates padding range to avoid filter-artifacts
cfg.artfctdef.eog.artpadding   = 0.4; % indicates range of data rejected beyond suprathresh.
cfg.artfctdef.eog.feedback     = 'yes';

cfgeye = ft_rejectartifact(cfg);
save(icfg.preprocfname,'cfgeye');

end

%%
function readdata(icfg)

load(icfg.preprocfname,'cfg');
cfg.channel         = {'MEG';'EEG057';'EEG058'};
cfg.bsfilter        = 'yes';
cfg.bsfreq          = [49.5,50.5];
cfg.padding         = 5;
cfg.blc             = 'yes';%'no';
cfg.detrend         = 'no';
data = ft_preprocessing(cfg);
save(icfg.FTdatafname,'data');
clear data

load(icfg.prepreyefname,'cfgeye');
cfg.trl            = cfgeye.artfctdef.eog.artifact;
cfg.trl(:,3)       = 0;
mycfg = cfg;
dateye              = ft_preprocessing(mycfg);
clear mycfg;

% rscfg.resamplefs    = 200;
% rscfg.detrend       = 'no';
% rscfg.blc           = 'yes';
% dateye              = ft_resampledata(rscfg,dateye);
save(icfg.FTeyedatafname,'dateye');
clear dateye

end

%% ALLOW REJECTION OF BAD COMPONENTS OR BAD CHANNELS (WHOLE TEMPORARY DATA SET)
function rawdatacomponentrejection(icfg)

close all;
display('################################');
display('#### eyemovement trials PCA ####');
display('#### ONLY reject VERY EXTREME outliers ####');
display('################################');
load(icfg.FTeyedatafname,'dateye');
eogchan = match_str(dateye.label,{'EEG057';'EEG058'});
megchan = match_str(dateye.label,ft_channelselection('MEG',dateye.label));

figeye = figure;
startover = 1;
while startover
    
    rcfg.method     = 'summary'; %'var';
    rcfg.channel    = 'MEG';
    rcfg.feedback   = 'yes';
    rcfg.keepchannel = 'yes';
    dateyedum          = ft_rejectvisual(rcfg,dateye);
    
    myvar = [];
    for ii=1:length(dateyedum.trial)
        myvar(ii,:) = nanvar(dateyedum.trial{ii}(megchan,:),[],2);
    end;
    figure;imagesc(myvar');ylabel('channels'),xlabel('trials');
    colorbar;title('Variance across trials and sensors before PCA');
    
    display('################################');
    display('#### running PCA decomposition of EYE trials ####');
    display('################################');
    
    cocfg.blc       = 'yes';
    cocfg.method    = 'pca';
    cocfg.channel   = {'MEG'};
    compeye         = ft_componentanalysis(cocfg,dateyedum);
    
    chandat = [];compdat = [];
    for iitrl=1:length(compeye.trial)
        chandat = cat(1,chandat,dateyedum.trial{iitrl}(eogchan,:)');
        compdat = cat(1,compdat,compeye.trial{iitrl}(1:4,:)');
    end;
    mycorr = corr([compdat,chandat]);
    
    display('PCA of raw temp data');
    
    figure(figeye);
    for pp = 1:8
        subplot(3,3,pp);
        if size(compeye.topo,1) == 271
            topoplot(compeye.topo(:,pp),'layout','CTF271Nott.lay','electrodes','off');
        elseif size(compeye.topo,1) == 270
            topoplot(compeye.topo(:,pp),'layout','CTF270Nott.lay','electrodes','off');
        elseif size(compeye.topo,1) == 274
            topoplot(compeye.topo(:,pp),'layout','CTFfil274.lay','electrodes','off');
        end;
    end;
    subplot(3,3,9);
    imagesc([1:6],[1:6],mycorr,[-1,1]);axis xy;
    set(gca,'XTick',[1:6],'YTick',[1:6]);
    set(gca,'XTickLabel',[compeye.label(1:4);dateye.label(eogchan)]);
    set(gca,'YTickLabel',[compeye.label(1:4);dateye.label(eogchan)]);
    colorbar;
    set(gcf,'Name','EOG detected');
    
    eyerccfg.component = input('Which bad PCA EYE components do you want to remove? E.g., 1 or [1,2] or [] ');
    eyerccfg.channel = {'MEG'};
    [datdum]   = ft_rejectcomponent(eyerccfg,compeye); %,datdum);
    
    % give feedback of effect of PCA removal
    myvar = [];
    for ii=1:length(datdum.trial)
        myvar(ii,:) = nanvar(datdum.trial{ii}(megchan,:),[],2);
    end;
    figure;imagesc(myvar');ylabel('channels'),xlabel('trials');
    colorbar;title('Variance across trials and sensors AFTER PCA');
    
    
    startover = input('Do you want to start all over again with trial rejection ? [0=no / 1=yes] ');
    clear cocfg datdum dateyedum;
    
end;

close all;
compeye = rmfield(compeye,'trial');
compeyeno = eyerccfg.component;
save(icfg.FTeyedatafname,'eyerccfg','compeye','compeyeno','dateye');
clear dateye

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% now work on the EPOCH DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
display('################################');
display('#### PCA on real data ####');
display('#### 1st ONLY reject VERY EXTREME outliers ####');
display('################################');

load(icfg.FTdatafname,'data');
eogchan = match_str(data.label,{'EEG057';'EEG058'});
megchan = match_str(data.label,ft_channelselection('MEG',data.label));

startover = 1;
while startover
    figmeg = figure;

    
    rcfg.method     = 'summary'; %'var';
    rcfg.channel    = 'MEG';
    rcfg.feedback   = 'yes';
    rcfg.keepchannel = 'yes';
    datdum           = ft_rejectvisual(rcfg,data);
    
    % give feedback of effect of PCA removal
    myvar = [];
    for ii=1:length(datdum.trial)
        myvar(ii,:) = nanvar(datdum.trial{ii}(megchan,:),[],2);
    end;
    figure;imagesc(myvar');ylabel('channels'),xlabel('trials');
    colorbar;title('Variance across trials and sensors BEFORE PCA');    
    
    cocfg.blc           = 'yes';
    cocfg.method        = 'pca';
    cocfg.channel       = {'MEG'};
    %cocfg.numcomponent  = length(data.label);
    comp                = ft_componentanalysis(cocfg,datdum);
    
    chandat = [];compdat = [];
    for iitrl=1:length(comp.trial)
        chandat = cat(1,chandat,datdum.trial{iitrl}(eogchan,:)');
        compdat = cat(1,compdat,comp.trial{iitrl}(1:4,:)');
    end;
    mycorr = corr([compdat,chandat]);
    
    display('#########################################');
    display('PCA of EPOCHED data');
    display('#########################################');
    figure(figmeg);
    for pp = 1:8
        subplot(3,3,pp);
        if size(comp.topo,1) == 271
            topoplot(comp.topo(:,pp),'layout','CTF271Nott.lay','electrodes','off');
        elseif size(comp.topo,1) == 270
            topoplot(comp.topo(:,pp),'layout','CTF270Nott.lay','electrodes','off');
        elseif size(comp.topo,1) == 274
            topoplot(comp.topo(:,pp),'layout','CTFfil274.lay','electrodes','off');
        end;
    end;
    subplot(3,3,9);
    imagesc([1:6],[1:6],mycorr,[-1,1]);axis xy;
    set(gca,'XTick',[1:6],'YTick',[1:6]);
    set(gca,'XTickLabel',[comp.label(1:4);data.label(eogchan)]);
    set(gca,'YTickLabel',[comp.label(1:4);data.label(eogchan)]);
    colorbar;
    set(gcf,'Name','all samples');
    
    noanswer = 1;
    while noanswer
        display('this last (real data) PCs (1),');
        whichcomponent = input('Do you want to use previous eye selective PCs (no, yes) ?  ','s');
        if strcmp(whichcomponent,'no')
            rccfg.component = input('Which bad PCA components do you want to remove? E.g., 1 or [1,2] or [] ');
            %rccfg.channel = {'MEG'};
            %[datdum]   = ft_rejectcomponent(rccfg,comp);
            noanswer = 0;
        elseif strcmp(whichcomponent,'yes') %use the eye components and apply to these data!
            %         [datdum]   = ft_rejectcomponent(eyerccfg,compeye,datdum);
            rccfg = eyerccfg;comp = compeye;
            noanswer = 0;
        else
            display('Give proper answer!!\n');
            noanswer = 1;
        end;
        
        mydat = datdum;
        mycomps = [1:size(comp.topo,1)];
        mycomps(rccfg.component) = [];
        for iitrl=1:length(datdum.trial)
            dummy = (datdum.trial{iitrl}(megchan,:)'*comp.topo);
            mydat.trial{iitrl} = (dummy(:,mycomps)*comp.unmixing(mycomps,:))';
            mydat.trial{iitrl}(eogchan,:) = datdum.trial{iitrl}(eogchan,:);
        end;
        datdum.trial = mydat.trial;
        
        % give feedback of effect of PCA removal
        myvar = [];
        for ii=1:length(datdum.trial)
            myvar(ii,:) = nanvar(datdum.trial{ii}(megchan,:),[],2);
        end;
        figure;imagesc(myvar');ylabel('channels'),xlabel('trials');
        colorbar;title('Variance across trials and sensors AFTER PCA');

        
    end;
    
    clc
    display('################################');
    display('#### rejectartefact on REAL DATA ####');
    display('#### This is your LAST CHANCE to remove artefacts!!!!!! ####');
    display('################################');
    rcfg.method     = 'summary'; %'var';
    rcfg.channel    = 'MEG';
    rcfg.feedback   = 'yes';
    rcfg.keepchannel = 'yes';
    datdum           = ft_rejectvisual(rcfg,datdum);
    
    startover = input('Do you want to start all over again with trial rejection ? [0=no / 1=yes] ');
    close all;
    
end;

data = datdum; clear datdum;
if isfield(comp,'trial'), comp = rmfield(comp,'trial'); end; %if using the compeye this wont exist
compno = rccfg.component;

save(icfg.FTdatafname,'data','rccfg','comp','compno');


end

function [directories,sessionfilenames,sessions] = collectsubjectinfo(dummy)

directories = {...
%     'S1'
%     'S2'
%     'S3';
%     'S4';
%     'S6';
%     'S7';
%     'S8';
%     'S9';
%     'S10';
%     'S11';
%     'S12';
%     'S13';
%     'S14';
%     'S16';
%     'S17';
%     'S18';
%     'S19';
%     'S20';
%     'S21';
%     'S24';
%     'S25';
%     'S22';
%     'S23';
    'S26';
    };

sessionfilenames = {...
%     '10987-010_MarkusBauer_20150717'
%     'S2-09467-008_MarkusBauer_20150721';
%     'S3-11133-01_MarkusBauer_20150724';
%     'S4-11139-001_MarkusBauer_20150721';
%     'mc09926-004_MarkusBauer_20150722';
%     'ia11147-001_MarkusBauer_20150723';
%     's8-10862-004_MarkusBauer_20150729';
%     'S9-11143-001_MarkusBauer_20150724';
%     'S10-11153-001_MarkusBauer_20150731';
%     '11151-001_MarkusBauer_20150728';
%     'S12-11126-003_MarkusBauer_20150729';
%     'S13-11024-007_MarkusBauer_20150804';
%     'S14-11013-001_MarkusBauer_20150804';
%     's16-11259-001_MarkusBauer_20151009';
%     's17-11245-001_MarkusBauer_20151004';
%     'S18-11261-001_MarkusBauer_20151009';
%     's19-11246-001_MarkusBauer_20151009';
%     '11254-001_MarkusBauer_20151007';
%     'jr11255-001_MarkusBauer_20151007';
%     's24-11260-001_MarkusBauer_20151009';
%     's25-11258-001_MarkusBauer_20151009';
%     's22-11270-001_MarkusBauer_20151012';
%     's23-11266-001_MarkusBauer_20151012';
    's26-11323-001_MarkusBauer_20151030';
    };

sessions = {...
%     1;
%     [2:5]; %tom
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     %[2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
%     [2:4];
    [2:4];
    };

end
