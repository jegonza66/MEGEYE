% (2018-03-26, Markus) Here is the battery for the SPM pipeline.
%% Set up path
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath(genpath('/home/juank/toolbox/fieldtrip-20170827'),'-end')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('/home/juank/toolbox/spm12/')
addpath('my_functions/')
addpath('markus/')
experimentdir   = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
rootmridir = [];%'d:\Projects\MSc_AttentionPrediction\MR_anatomicals\';
directories     = {'/13229-001'};
sessionfilenames= {'13229-001_MarkusBauer_20180319'};
sessions        = {[4:6]};

load('labels_CTF269Nott','labels');

subjectinfo.directories         = directories;
subjectinfo.sessionfilenames    = sessionfilenames;

%% First steps of analysis (PREPARE DATA FOR COVARIANCE)
% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess))];
        cfg.dataset = [cfg.directory,cfg.myfname,'.ds'];
    
        if exist(cfg.dataset,'dir')
            cfg.datafname       = [cfg.directory,'/tmpcovariance/data_',        cfg.myfname,'.mat'];
            cfg.datafname_weyech= [cfg.directory,'/tmpcovariance/data_weyech_', cfg.myfname,'.mat'];
            cfg.preprocfname    = [cfg.directory,'/tmpcovariance/cfg_',         cfg.myfname,'.mat'];

            cfg.eyemap_datafname       = [cfg.directory,'/tmpcovariance/eyemap_data_', cfg.myfname,'.mat'];
            cfg.eyemap_datafname_weyech= [cfg.directory,'/tmpcovariance/eyemap_data_weyech_', cfg.myfname,'.mat'];
            cfg.eyemap_preprocfname    = [cfg.directory,'/tmpcovariance/eyemap_cfg_',  cfg.myfname,'.mat'];

            if ~exist([cfg.directory,'/tmpcovariance'],'dir')
                mkdir([cfg.directory,'/tmpcovariance']);
            end;

            % DO THE EPOCHING & BLINK-DETECTION
            % Setup eventmarks (saves a configuration file).
            fun_epoching(cfg);

            % Read data and perform some preprocessing (saves a data file).
            fun_readdata_simple(cfg);
        end;
    end;

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
%     cfg.sequence            = [cfg.directory,'/tmpcovariance/sequence_', cfg.myfname,'.mat'];
%     cfg.eyemap_sequence     = [cfg.directory,'/tmpcovariance/eyemap_sequence_', cfg.myfname,'.mat'];
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
    cfg.datafname               = [cfg.directory,'/tmpcovariance/data_', cfg.myfname,'.mat'];
    cfg.datafname_weyech        = [cfg.directory,'/tmpcovariance/data_weyech_', cfg.myfname,'.mat'];
    cfg.preprocfname            = [cfg.directory,'/tmpcovariance/cfg_',  cfg.myfname,'.mat'];

    cfg.eyemap_datafname        = [cfg.directory,'/tmpcovariance/eyemap_data_', cfg.myfname,'.mat'];
    cfg.eyemap_datafname_weyech = [cfg.directory,'/tmpcovariance/eyemap_data_weyech_', cfg.myfname,'.mat'];
    cfg.eyemap_preprocfname     = [cfg.directory,'/tmpcovariance/eyemap_cfg_',  cfg.myfname,'.mat'];

    cfg.sequence            = [cfg.directory,'/tmpcovariance/sequence_', cfg.myfname,'.mat'];
    cfg.eyemap_sequence     = [cfg.directory,'/tmpcovariance/eyemap_sequence_', cfg.myfname,'.mat'];
    
    fun_replace_eventmarkers(cfg.eyemap_datafname,          cfg.eyemap_datafname,       cfg.eyemap_sequence)
    fun_replace_eventmarkers(cfg.eyemap_datafname_weyech,   cfg.eyemap_datafname_weyech,cfg.eyemap_sequence)
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

        cfg.myfname = ['eyemap_data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];
        fun_rawdatacomponentrejection_v3(datafname,datafname,'previous_unmixing',{comp,compno});
    end
    
%% Calculate covariance matrix and noise covariance matrix
% supertimelock_calccovar.m

% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];
        
        load(datafname,'data');
        
        load('labels_CTF269Nott','labels');
        scfg.channel = labels;
        data = ft_selectdata(scfg,data);
        
        % data.grad
        %      balance: [1x1 struct]
        %      chanori: [296x3 double]
        %      chanpos: [296x3 double]
        %     chantype: {296x1 cell}
        %     chanunit: {296x1 cell}
        %      coilori: [583x3 double]
        %      coilpos: [583x3 double]
        %     coordsys: 'ctf'
        %        label: {296x1 cell}
        %          tra: [296x583 double]
        %         type: 'ctf275'
        %         unit: 'cm'        
        
        % Shouldn't be 269???
        if isfield(data,'label') && (length(data.label) > 274)
            data.label = data.label(1:274);
        end;
        
        timelockfile = [cfg.directory,'tmpcovariance',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
        tmlkcfg             = [];
        tmlkcfg.vartrllength= 0;
        tmlkcfg.channel     = 'MEG';
        tmlkcfg.lpfilter    = 'yes';
        tmlkcfg.lpfreq      = 45;
        tmlkcfg.keeptrials  = 'yes';
        tmlkcfg.covariance  = 'yes';
        
        timelock = ft_timelockanalysis(tmlkcfg,data); 
        timelock.grad       = ft_convert_units(timelock.grad,'m');
        save(timelockfile,'timelock');
        
        timelockfile = [cfg.directory,'tmpcovariance',filesep,'tmlklcmvhifrqfilt_',num2str(isess),'.mat'];
        tmlkcfg             = [];
        tmlkcfg.vartrllength= 0;
        tmlkcfg.channel     = 'MEG';
        tmlkcfg.hpfilter    = 'yes';
        tmlkcfg.hpfreq      = 40;
        tmlkcfg.keeptrials  = 'yes';
        tmlkcfg.covariance  = 'yes';

        timelock = ft_timelockanalysis(tmlkcfg,data); 
        timelock.grad       = ft_convert_units(timelock.grad,'m');
        save(timelockfile,'timelock');
    end;
    
% end;

%% And the inverse solution...
% super_calculate_LCMVfilter_SAM.m
% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
% isess=1;
        headmodelfile = [cfg.directory,'tmpcovariance',filesep,'fwdmdlspm12raw_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        load(headmodelfile);
        positions = forward.forward.mesh.vert;

        timelockfile = [cfg.directory,'tmpcovariance',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
        load(timelockfile);
        
%         [~,unselchans] = match_str({'MLO12';'MLO23';'MRO41';'MRO42'},timelock.label);
%         if ~isempty(unselchans)
%             if (length(size(timelock.cov))==2)
%                     timelock.cov(unselchans,:) = [];
%                     timelock.cov(:,unselchans) = [];
%                     timelock.avg(unselchans,:) = [];
%                     timelock.dof(unselchans,:) = [];
%                     timelock.var(unselchans,:) = [];
%                     labels = timelock.label(setdiff([1:length(timelock.label)]',unselchans));
%                     timelock.label = labels;
%             elseif (length(size(timelock.cov))==3)
%                     timelock.cov(:,unselchans,:) = [];
%                     timelock.cov(:,:,unselchans) = [];
%                     timelock.avg(unselchans,:) = [];
%                     timelock.dof(unselchans,:) = [];
%                     timelock.var(unselchans,:) = [];
%                     labels = timelock.label(setdiff([1:length(timelock.label)]',unselchans));
%                     timelock.label = labels;
%             end            
%         end;
        
        gsrcfg.vol          = forward.forward.vol;
        gsrcfg.grid.pos     = positions;
        gsrcfg.feedback     = 'none';
        gsrcfg.channel      = {'MEG'}; %;'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
        gsrcfg.reducerank   = 2;
        gsrcfg.sourceunits  = 'm';
        gsrcfg.grad         = ft_convert_units(timelock.grad,'m');

        leadfieldfile = [cfg.directory,'tmpcovariance',filesep,'leadfieldBNDEM_',num2str(sessions{isubject}(isess)),'.mat'];
        if ~exist(leadfieldfile,'file')
            gsrcfg.inwardshift = 0;
            grid = ft_prepare_leadfield(gsrcfg,timelock);
            save(leadfieldfile,'grid');
        else
            load(leadfieldfile);
        end;
% And also it says: 
% 0 dipoles inside, 8196 dipoles outside brain
% Is it ok? Is it just over the brain (i.e. the cortex) or is it further
% away (i.e. AN ERROR)?
% SOLVED! I only needed to change the units of everything to m.

% figure(1);clf
%     plot3(gsrcfg.vol.bnd.pos(:,1),gsrcfg.vol.bnd.pos(:,2),gsrcfg.vol.bnd.pos(:,3),'.')
%     hold on;plot3(gsrcfg.grid.pos(:,1),gsrcfg.grid.pos(:,2),gsrcfg.grid.pos(:,3),'.k')
%     hold on;plot3(gsrcfg.grad.chanpos(:,1),gsrcfg.grad.chanpos(:,2),gsrcfg.grad.chanpos(:,3),'.r')

% UNFILTERED LCMV BEAMFORMER
    srcfg.method        = 'sam';
    srcfg.feedback      = 'none';
    srcfg.channel       = {'MEG'}; %'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
    srcfg.vol           = gsrcfg.vol;
%         srcfg.meansphereorigin = nanmean(srcfg.vol.bnd.pnt,1);
    srcfg.meansphereorigin = nanmean(srcfg.vol.bnd.pos,1);
                        % Reference to non-existent field 'pnt'.
                        % It should be pos, because it makes sense with
                        % the name (tri comes from the triangles), and
                        % because .pnt often comes along with the .tri
                        % in the manuals and blogs.
    srcfg.grad          = gsrcfg.grad;
    srcfg.parallel      = 'no';
    srcfg.keepfilter    = 'yes';
    srcfg.fixedori      = 'spinning';
    srcfg.grid          = grid;
%     srcfg.keeptrials    = 'yes';
    
    figure(1); clf
        plot3(srcfg.grid.pos(:,1),srcfg.grid.pos(:,2),srcfg.grid.pos(:,3),'b.')
        hold on; plot3(srcfg.vol.bnd.pos(:,1),srcfg.vol.bnd.pos(:,2),srcfg.vol.bnd.pos(:,3),'r.')
        hold on; plot3(srcfg.grad.chanpos(:,1),srcfg.grad.chanpos(:,2),srcfg.grad.chanpos(:,3),'ko')
%     keyboard
    
    if (length(size(timelock.cov))==2)
        srcfg.lambda        = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
    elseif (length(size(timelock.cov))==3)
        srcfg.lambda        = nanmean(abs(diag(squeeze(mean(timelock.cov)))))/50; %nanmean(abs(timelock.cov(:)));
%         srcfg.lambda = nan(1,size(timelock.cov,1));
%         for i=1:size(timelock.cov,1)
%             srcfg.lambda(i)        = nanmean(abs(diag(squeeze(timelock.cov(i,:,:)))))/50; %nanmean(abs(timelock.cov(:)));
%         end
    end
    source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));
    %filter = source.avg.filter;

    filterfile = [cfg.directory,'tmpcovariance',filesep,'LCMVBNDEMnofiltSAMhires_',num2str(sessions{isubject}(isess))];
    save(filterfile,'source');
    clear source timelock;

% FILTER CALCULATION - TIMELOCK - high-pass filtered
    timelockfile = [cfg.directory,'tmpcovariance',filesep,'tmlklcmvhifrqfilt_',num2str(isess),'.mat'];
    load(timelockfile);

%         [~,unselchans] = match_str({'MLO12';'MLO23';'MRO41';'MRO42'},timelock.label);
%         if ~isempty(unselchans)
%             timelock.cov(unselchans,:) = [];
%             timelock.cov(:,unselchans) = [];
%             timelock.avg(unselchans,:) = [];
%             timelock.dof(unselchans,:) = [];
%             timelock.var(unselchans,:) = [];
%             labels = timelock.label(setdiff([1:length(timelock.label)]',unselchans));
%             timelock.label = labels;
%         end;        

    figure(1); clf
        plot3(srcfg.grid.pos(:,1),srcfg.grid.pos(:,2),srcfg.grid.pos(:,3),'b.')
        hold on; plot3(srcfg.vol.bnd.pos(:,1),srcfg.vol.bnd.pos(:,2),srcfg.vol.bnd.pos(:,3),'r.')
        hold on; plot3(srcfg.grad.chanpos(:,1),srcfg.grad.chanpos(:,2),srcfg.grad.chanpos(:,3),'ko')
%     keyyboard
    
    if (length(size(timelock.cov))==2)
        srcfg.lambda        = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
    elseif (length(size(timelock.cov))==3)
        srcfg.lambda        = nanmean(abs(diag(squeeze(mean(timelock.cov)))))/50; %nanmean(abs(timelock.cov(:)));
%         srcfg.lambda = nan(1,size(timelock.cov,1));
%         for i=1:size(timelock.cov,1)
%             srcfg.lambda(i)        = nanmean(abs(diag(squeeze(timelock.cov(i,:,:)))))/50; %nanmean(abs(timelock.cov(:)));
%         end
    end
    source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));
    %filter = source.avg.filter;

    filterfile = [cfg.directory,'tmpcovariance',filesep,'LCMVBNDEMhifrqfiltSAMhires_',num2str(sessions{isubject}(isess))];
    save(filterfile,'source');
    clear source timelock grid srcfg gsrcfg;
        
        
    end;
    
% end;

%% Project into ROIs estimating the PCA within each ROI.
% super_freqanalysis_source_ROISAM_v2hires.m
% cd(experimentdir)
% mysurfhi = load([experimentdir,'analysisroutines\canonical_surf_MNI.mat']);
% mysurflo = load([experimentdir,'analysisroutines\canonical_surf_MNI_413.mat']);
mysurfhi = load('markus/canonical_surf_MNI.mat');
mysurflo = load('markus/canonical_surf_MNI_413.mat');

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);

% isess = 1;
for isess = 1:length(sessions{isubject})
    cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
    datafname = [cfg.directory,'tmpcovariance',filesep,cfg.myfname];

    load(datafname,'data');
    
    load('labels_CTF269Nott','labels');
    scfg.channel = labels;
    data = ft_selectdata(scfg,data);

%     data.hdr.grad   = ft_convert_units(data.hdr.grad,'m');
%     data.grad       = ft_convert_units(data.grad,'m');
    figure(1); clf; plot3(data.grad.chanpos(:,1),data.grad.chanpos(:,2),data.grad.chanpos(:,3),'b.')
    
% LOW frequencies
    filterfile = [cfg.directory,'tmpcovariance',filesep,'LCMVBNDEMnofiltSAMhires_',num2str(sessions{isubject}(isess))];
    load(filterfile);
    hold on; plot3(source.pos(:,1),source.pos(:,2),source.pos(:,3),'r.')
    % Sources are in meters
    
    searchregion = 7; %7mm  % Be careful!!! Everything is in meters now!!!
                            % Not everything, the surf.XXX are in cm
                            
    figure(1); clf;     plot3(mysurfhi.surf.pnt(:,1),mysurfhi.surf.pnt(:,2),mysurfhi.surf.pnt(:,3),'b.')
            hold on;    plot3(mysurflo.surf.pnt(:,1),mysurflo.surf.pnt(:,2),mysurflo.surf.pnt(:,3),'ro')
            
    sconfig = struct();
    sconfig.filter              = findbeamfcoeffs4ROIlowreshiresgrid(...
                            mysurflo.surf,mysurfhi.surf,source,searchregion);
    % source.avg.filter is cell array of Nsources (8196) vectors with dim = Nsensors (269),
    % i.e. it has the projection of the sensors into each source
    % sconfig.filter is cell array of Nsources arrays with dim = 
    % Nneighbor sources (~12) x Nsensors,
    sconfig.sourcelabel         = source.cfg.previous.channel;
    sconfig.preproc.lpfilter    = 'yes'; 
    sconfig.preproc.lpfreq      = 40;        
    clear filter source

    [datdum] = projectdatathrough_allfilters_ROIPCAlatest(data,sconfig);

    dattmp = data; data = datdum;
    save([cfg.directory,'sources',filesep,...
        'data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'_src413nofilt.mat'],'data');
    data = dattmp; clear dattmp
    
%     cfg.freqlowfname = [cfg.directory,'sources',filesep,'freqLowhiressess',num2str(isess)];
%     freqanalysisLOW_predatt_source(cfg,datdum);
%     clear datdum;

% HIGH frequencies
    filterfile = [cfg.directory,'tmpcovariance',filesep,'LCMVBNDEMhifrqfiltSAMhires_',num2str(sessions{isubject}(isess))];
    load(filterfile);

    searchregion = 7; %mm %% Be careful!!! Everything is in meters now!!!
    sconfig = struct();
    sconfig.filter          = findbeamfcoeffs4ROIlowreshiresgrid(...
                        mysurflo.surf,mysurfhi.surf,source,searchregion);
    sconfig.sourcelabel     = source.cfg.previous.channel;
    sconfig.preproc.hpfilter= 'yes'; 
    sconfig.preproc.hpfreq  = 35;
    clear filter source

    [datdum] = projectdatathrough_allfilters_ROIPCAlatest(data,sconfig);

    dattmp = data; data = datdum;
    save([cfg.directory,'sources',filesep,...
            'data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'_src413freqHigh.mat'],'data');
    data = dattmp; clear dattmp
    
%     cfg.freqhighfname = [cfg.directory,'sources',filesep,'freqHighhiressess',num2str(isess)];
%     freqanalysisHIGH_predatt_source(cfg,datdum);
%     clear source datdum data;

end;
% end;

%% Before following the analysis I'm need to join the three parts of the experiment.
% I plan to do it straightforward, as I did it before, but we shall see
% what happen with the sources.

% super_freqanalysis_source_ROISAM_v2hires.m
% cd(experimentdir)
% mysurfhi = load([experimentdir,'analysisroutines\canonical_surf_MNI.mat']);
% mysurflo = load([experimentdir,'analysisroutines\canonical_surf_MNI_413.mat']);
mysurfhi = load('markus/canonical_surf_MNI.mat');
mysurflo = load('markus/canonical_surf_MNI_413.mat');

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);

filetype = {'src413nofilt'};
for isess = 1:length(sessions{isubject})
%     cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'_src413nofilt.mat'];
    cfg.myfname = [filetype{1},num2str(isess)];
    datafname = [cfg.directory,'sources',filesep,cfg.myfname];

    load(datafname,'data');
    data2merge{isess} = data;
end
data        = ft_appenddata([], data2merge{1}, data2merge{2}, data2merge{3});
cfg.myfname = ['data_',sessionfilenames{isubject},'_',filetype{1},'_merged.mat'];
save([cfg.directory,'/sources/',cfg.myfname],'data')

%% For the synchronization I need the eye channels, 
% so I'm going to merge here 
isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

display('________________________________________________');
display(['subject: ',num2str(isubject)]);

filetypes = {'data_weyech_'};
for isess = 1:length(sessions{isubject})
    cfg.myfname = [filetypes{1},sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];

    load([cfg.directory,'/tmpcovariance/',cfg.myfname]);
    data2merge{isess} = data;
end
data = ft_appenddata([], data2merge{1}, data2merge{2}, data2merge{3});
save([cfg.directory,'/sources/',cfg.myfname(1:end-7),'_merged.mat'],'data')

%% Synchronize the MEG time with the ET time
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
addpath('/home/juank/toolbox/fieldtrip-20170827')
ft_defaults;
addpath('my_functions/')
experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319'};

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

% Load merged files
cfg.myfname = ['data_weyech_',sessionfilenames{isubject},'_merged.mat'];
load([cfg.directory,'sources/',cfg.myfname],'data')

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

data.fsample = round(mean(1./diff(data.time{1})));
[eyedata,xc_xy] = fun_estimating_lag(data,eyedata,0);
fprintf('Total fixations = %d\n',sum([eyedata.Nfix]))

% Load task data
cfg.myfname = ['data_',sessionfilenames{isubject},'_src413nofilt_merged.mat'];
load([cfg.directory,'sources/',cfg.myfname],'data')
data.fsample = round(mean(1./diff(data.time{1})));

recfg = [];
recfg.filterNonStimFix = 1;
recfg.epochLimits_ms= [-200 250]; % In ms,
recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[dataEpoched,recfg] = fun_reepoch_data(recfg,data,eyedata);

cfg.myfname = ['dataEpoched_',sessionfilenames{isubject},'_src413nofilt_merged.mat'];
save([cfg.directory,'sources/',cfg.myfname],'dataEpoched','recfg')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Instead of projecting I want to keep all the points just to see how it looks



