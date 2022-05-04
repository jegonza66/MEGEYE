
function super_preprocessPCAFT(mysubjects)

%[directories,sessionfilenames] = collectsubjectinfo;
addpath d:\matlab\donders\Ole\fieldtrip\

experimentdir = 'd:\Projects\MichaelSchmid\';
directories = {%'Data\';
               %'Data\';
               %'Data\';
               %'Data\';
               'Data\';
               'Data\';
               'Data\';
               };
sessionfilenames = {%'thetatest_MarkusBauer_20180206_03.ds'; %Samy stimuli close together, small and RF compatible
                    %'thetatest_MarkusBauer_20180206_04.ds'; %Samy stimuli further apart, small and RF compatible
                    %'thetatest_MarkusBauer_20180206_01.ds'; %Michael stimuli "close" together, large and relatively further apart
                    %'thetatest_MarkusBauer_20180206_02.ds'; %Michael stimuli further apart, large and relatively further apart
                    'thetatest_MarkusBauer_20180205_01.ds';
                    'thetatest_MarkusBauer_20180205_02.ds';
                    'thetatest_MarkusBauer_20180205_03.ds';
                    };

%[directories,sessionfilenames] = collectsubjectinfo;
% if nargin<1
%     mysubjects = [1:length(directories)];
% end;
% mysubjects

for j = 1:length(sessionfilenames)
    
    cfg.directory = [experimentdir,directories{j}];
    
    display('________________________________________________');
    display(['subject: ',num2str(j)]);
    cfg.myfname = sessionfilenames{j}([1:(end-3)]);
    
    cfg.dataset = [cfg.directory,filesep,sessionfilenames{j}];
    
    if exist(cfg.dataset,'dir')
        
        cfg.datafname = [cfg.directory,'\matfiles\data_',cfg.myfname,'.mat'];
        cfg.preprocfname = [cfg.directory,'\matfiles\cfg_',cfg.myfname,'.mat'];
        if ~exist([cfg.directory,'\matfiles'],'dir')
            mkdir([cfg.directory,'\matfiles']);
        end;
        
        epoching(cfg);
        readdata(cfg);
        rawdatacomponentrejection(cfg);
        
    end;
    
end;

end



%% DO THE EPOCHING & BLINK-DETECTION

function epoching(icfg)


cfg.dataset                    = icfg.dataset;
cfg.continuous                 = 'yes';

%% define trials / epochs
% cfg.trialdef.eventvalue      = []; %name of the triggers in brain vision format
% cfg.trialdef.eventtype          = {'Relevant';'Irrelevant'};
cfg.trialdef.eventtype          = 'gui';
cfg.trialdef.prestim            = 2;
cfg.trialdef.poststim           = 4;
cfg = ft_definetrial(cfg);

%cfg.trl([1,end],:) = [];
save(icfg.preprocfname,'cfg');

end

function readdata(icfg)

load(icfg.preprocfname,'cfg');
cfg.channel         = {'meg';'megref'};
%cfg.channel         = {'MEG'};
cfg.bsfilter        = 'yes';
cfg.bsfreq          = [49.5,50.5];
cfg.padding         = 7.168;
cfg.blc             = 'yes';%'no';
cfg.detrend         = 'no';
data = ft_preprocessing(cfg);

ccfg.resamplefs = 600;ccfg.demean = 'yes';ccfg.detrend = 'yes';
[data] = ft_resampledata(ccfg, data);

%now the data with reference sensors (aligned with the MEG sensors in time)
%and ocular artefact removed are transformed to, e.g. 1st grad
gcfg.gradient = 'G3BR'; %'none'
data = ft_denoise_synthetic(gcfg,data);

load('labels_CTF269Nott','labels');
scfg.channel = labels;data = ft_selectdata(scfg,data);
save(icfg.datafname,'data');
 
% %% specify preprocessing options
% cfg.channel                     = {'MEG'};
% % filtering - should be done before epoching ideally
% % line noise - bandstop filter
% cfg.bsfilter                    = 'yes';
% cfg.bsfreq                      = [49,51];
% cfg.padding                     = 10;
% cfg.blc                         = 'yes';%'no';
% cfg.detrend                     = 'no';
% cfg.artfctdef.feedback          = 'yes';
% cfg.artfctdef.reject            = 'partial';
% cfg.artfctdef.minaccepttim      = 0.1;
% cfg.artfctdef.type              = {'eog'};
% 
% cfg.artfctdef.eog.bpfilter     = 'yes';
% cfg.artfctdef.eog.bpfilttype   = 'but';
% cfg.artfctdef.eog.bpfreq       = [1 15];
% cfg.artfctdef.eog.bpfiltord    = 4;
% cfg.artfctdef.eog.hilbert      = 'yes';
% cfg.artfctdef.eog.channel      = {'EOG'};
% cfg.artfctdef.eog.cutoff       = 3;   % z-value at which to threshold
% cfg.artfctdef.eog.trlpadding   = 0;   % indicates range to search 4 artifact
% cfg.artfctdef.eog.fltpadding   = 1;   % indicates padding range to avoid filter-artifacts
% cfg.artfctdef.eog.artpadding   = 0.4; % indicates range of data rejected beyond suprathresh.
% cfg.artfctdef.eog.feedback     = 'yes';
% 
% cfgeye = ft_rejectartifact(cfg);
% fnamestring = [datapath,'\cfgEYEsess'];
% save(fnamestring,'cfgeye');
% close all;

end


%% ALLOW REJECTION OF BAD COMPONENTS OR BAD CHANNELS (WHOLE TEMPORARY DATA SET)

function rawdatacomponentrejection(icfg)

load(icfg.datafname,'data');

% scfg.channel = {'all', '-EOG'};
% data = ft_selectdata_new(scfg,data);

startover = 1;
while startover
    
    rcfg.method     = 'summary';
    rcfg.channel    = 'all';
    rcfg.keepchannel = 'nan';
    rcfg.feedback   = 'yes';
    datdum           = ft_rejectvisual(rcfg,data);   

    cocfg.blc           = 'yes';
    cocfg.method        = 'pca';
    cocfg.channel       = {'all'};
    cocfg.numcomponent  = length(data.label);
    comp                = ft_componentanalysis(cocfg,datdum);
    
    display('PCA of raw temp data');
    figure;
    for pp = 1:8
        subplot(2,4,pp);
        if size(comp.topo,1) == 269
            topoplot(comp.topo(:,pp),'layout','CTF269Nott.lay','electrodes','off');
        elseif size(comp.topo,1) == 270
            topoplot(comp.topo(:,pp),'layout','CTF270Nott.lay','electrodes','off');
        end;
    end;   
    
        rccfg.component = input('Which bad PCA components do you want to remove? E.g., 1 or [1,2] or [] ');
        [datdum]   = ft_rejectcomponent(rccfg,comp,datdum);        
    
    rcfg.method     = 'summary';
    rcfg.channel    = 'all';
    rcfg.keepchannel = 'nan';
    rcfg.feedback   = 'yes';
    datdum           = ft_rejectvisual(rcfg,datdum);   
    
    startover = input('Do you want to start all over again with trial rejection ? [0=no / 1=yes] ');
    clear cocfg;
    close all;
    
end;

data = datdum; clear datdum;
comp = rmfield(comp,'trial'); 
compno = rccfg.component;

save(icfg.datafname,'data','rccfg','comp','compno');

end


% function [directories,sessionfilenames] = collectsubjectinfo(dummy)
% 
% directories = {...
%     '1001'
%     '1002'
%     '1003'
%     '1005'
%     '1007'
%     '1009'
%     '1010'
%     '1011'
%     '1014'
%     '1016'
%     '1021'
%     '1024'
%     '1032'
%     '1033'
%     '1041'
%     '1047'
%     '1048'
%     '1054'
%     '1061'
%     '1068'
%     '1072'
%     '1080'
% %     '1081'
% %     '1088'
% %     '1089'
% %     '1091'
% %     '1093'
% %     '1094'
%     };
% 
% sessionfilenames = {...
%     '1001_relmod.ds'
%     '1002_relmod.ds'
%     '1003_relmod.ds'
%     '1005_relmod.ds'
%     '1007_relmod.ds'
%     '1009_relmod.ds'
%     '1010_relmod.ds'
%     '1011_relmod.ds'
%     '1014_relmod.ds'
%     '1016_relmod.ds'
%     '1021_relmod.ds'
%     '1024_relmod.ds'
%     '1032_relmod.ds'
%     '1033_relmod.ds'
%     '1041_relmod.ds'
%     '1047_relmod.ds'
%     '1048_relmod.ds'
%     '1054_relmod.ds'
%     '1061_relmod.ds'
%     '1068_relmod.ds'
%     '1072_relmod.ds'
%     '1080_relmod.ds'
% %     '1081_relmod.ds'
% %     '1088_relmod.ds'
% %     '1089_relmod.ds'
% %     '1091_relmod.ds'
% %     '1093_relmod.ds'
% %     '1094_relmod.ds'
%     };
% 
% end
% 
