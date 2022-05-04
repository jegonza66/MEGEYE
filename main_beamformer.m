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
experimentdir   = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
rootmridir = [];%'d:\Projects\MSc_AttentionPrediction\MR_anatomicals\';
directories     = {'/13229-001'};
sessionfilenames= {'13229-001_MarkusBauer_20180319'};
sessions        = {[4:6]};

load('labels_CTF269Nott','labels');

subjectinfo.directories         = directories;
subjectinfo.sessionfilenames    = sessionfilenames;

%% Calculate covariance matrix and noise covariance matrix
% supertimelock_calccovar.m

% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = ['data_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        datafname = [cfg.directory,'matfiles',filesep,cfg.myfname];
        
        load(datafname,'data');
        
        % This data
        % First, I want to preprocess a little bit the data, at least
        % removing the ICA components.
        
        % Shouldn't be 269???
        if isfield(data,'label') && (length(data.label) > 274)
            data.label = data.label(1:274);
        end;
        
        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
        tmlkcfg             = [];
        tmlkcfg.vartrllength= 0;
        tmlkcfg.channel     = 'MEG';
        tmlkcfg.lpfilter    = 'yes';
        tmlkcfg.lpfreq      = 45;
        tmlkcfg.keeptrials  = 'no';
        tmlkcfg.covariance  = 'yes';
        
        timelock = ft_timelockanalysis(tmlkcfg,data); 
        save(timelockfile,'timelock');
        
        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvhifrqfilt_',num2str(isess),'.mat'];
        tmlkcfg             = [];
        tmlkcfg.vartrllength= 0;
        tmlkcfg.channel     = 'MEG';
        tmlkcfg.hpfilter    = 'yes';
        tmlkcfg.hpfreq      = 40;
        tmlkcfg.keeptrials  = 'no';
        tmlkcfg.covariance  = 'yes';
        
        timelock = ft_timelockanalysis(tmlkcfg,data); 
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
    
%     for isess = 1:length(sessions{isubject})
isess=1;
        headmodelfile = [cfg.directory,'matfiles',filesep,'fwdmdlspm12raw_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
        load(headmodelfile);
        positions = forward.forward.mesh.vert;

        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
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
        
        gsrcfg.vol          = forward.forward.vol;
        gsrcfg.grid.pos     = positions;
        gsrcfg.feedback     = 'none';
        gsrcfg.channel      = {'MEG'}; %;'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
        gsrcfg.reducerank   = 2;
        gsrcfg.sourceunits  = 'mm';
        gsrcfg.grad         = timelock.grad;
        
        leadfieldfile = [cfg.directory,'matfiles',filesep,'leadfieldBNDEM_',num2str(sessions{isubject}(isess)),'.mat'];
        if ~exist(leadfieldfile,'file')
            gsrcfg.inwardshift = 0;
            grid = ft_prepare_leadfield(gsrcfg,timelock);
            save(leadfieldfile,'grid');
        else
            load(leadfieldfile);
        end;
% Check these WARNINGs, Is it just a change of name? Or is it doing
% something else?
% Warning: use cfg.headmodel instead of cfg.vol 
% Warning: use cfg.unit instead of cfg.sourceunits 
% Warning: The field cfg.unit is deprecated, pleae use cfg.grid.unit        

% And also it says: 
% 0 dipoles inside, 8196 dipoles outside brain
% Is it ok? Is it just over the brain (i.e. the cortex) or is it further
% away (i.e. AN ERROR)?


% hold on;plot3(gsrcfg.grid.pos(:,1),gsrcfg.grid.pos(:,2),gsrcfg.grid.pos(:,3),'.k')
% hold on;plot3(gsrcfg.grad.pos(:,1),gsrcfg.grad.pos(:,2),gsrcfg.grad.pos(:,3),'.r')
% hold on;plot3(gsrcfg.grad.chanpos(:,1),gsrcfg.grad.chanpos(:,2),gsrcfg.grad.chanpos(:,3),'.r')

% gridd = ft_convert_units(gsrcfg.grid,'cm')
        %% UNFILTERED LCMV BEAMFORMER
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
        srcfg.lambda        = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
        
        source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));
        %filter = source.avg.filter;
        
        filterfile = [cfg.directory,'matfiles',filesep,'LCMVBNDEMnofiltSAMhires_',num2str(sessions{isubject}(isess))];
        save(filterfile,'source');
        clear source timelock;
        
% I get again this warning!!
% Warning: use cfg.headmodel instead of cfg.vol        
% number of dipoles inside  brain: 0
% number of dipoles outside brain: 8196

        %% FILTER CALCULATION - TIMELOCK - high-pass filtered
        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvhifrqfilt_',num2str(isess),'.mat'];
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
        
        srcfg.lambda = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
        source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));
        %filter = source.avg.filter;
        
        filterfile = [cfg.directory,'matfiles',filesep,'LCMVBNDEMhifrqfiltSAMhires_',num2str(sessions{isubject}(isess))];
        save(filterfile,'source');
        clear source timelock grid srcfg gsrcfg;
        
        
%     end;
    
% end;