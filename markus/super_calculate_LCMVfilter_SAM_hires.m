
function super_calculate_LCMVfilter_SAM(mysubjects)


switchtofieldtrip20150923
experimentdir = 'd:\Projects\MSc_AttentionPrediction\MEG\';
[directories,sessionfilenames,sessions] = collectsubjectinfo;

if nargin<1
    mysubjects = [1:length(directories)];
end;
mysubjects

for isubject = mysubjects
    
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        
        headmodelfile = [cfg.directory,'matfiles',filesep,'fwdmdlspm8raw_',sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.mat'];
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
        
        gsrcfg.vol = forward.forward.vol;
        gsrcfg.grid.pos = positions;
        gsrcfg.feedback    = 'none';
        gsrcfg.channel     = {'MEG'}; %;'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
        gsrcfg.reducerank  = 2;
        gsrcfg.sourceunits = 'mm';
        gsrcfg.grad = timelock.grad;
        
        leadfieldfile = [cfg.directory,'matfiles',filesep,'leadfieldBNDEM_',num2str(sessions{isubject}(isess)),'.mat'];
        if ~exist(leadfieldfile,'file')
            gsrcfg.inwardshift = 0;
            grid = ft_prepare_leadfield(gsrcfg,timelock);
            save(leadfieldfile,'grid');
        else
            load(leadfieldfile);
        end;
        
        
        %% UNFILTERED LCMV BEAMFORMER
        srcfg.method        = 'sam';
        srcfg.feedback      = 'none';
        srcfg.channel       = {'MEG'}; %'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
        srcfg.vol           = gsrcfg.vol;
        srcfg.meansphereorigin = nanmean(srcfg.vol.bnd.pnt,1);
        srcfg.grad          = gsrcfg.grad;
        srcfg.parallel      = 'no';
        srcfg.keepfilter    = 'yes';
        srcfg.fixedori      = 'spinning';
        srcfg.grid          = grid;
        srcfg.lambda = nanmean(abs(diag(timelock.cov)))/50; %nanmean(abs(timelock.cov(:)));
        
        source = ft_source2sparse(ft_sourceanalysis(srcfg,timelock));
        %filter = source.avg.filter;
        
        filterfile = [cfg.directory,'matfiles',filesep,'LCMVBNDEMnofiltSAMhires_',num2str(sessions{isubject}(isess))];
        save(filterfile,'source');
        clear source timelock;
        
        
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
        
        
    end;
    
end;

end


function [directories,sessionfilenames,sessions] = collectsubjectinfo(dummy)

directories = {...
    'S1';
    'S2';
    'S3';
    'S4';
    'S6';
    'S7';
    'S8';
    'S9';
    'S10';
    'S11';
    'S12';
    'S13';
    'S14';
    'S16';
    'S17';
    'S18';
    'S19';
    'S20';
    'S21';
    'S22';
    'S23';
    'S24';
    'S25';
    'S26';
    };

sessionfilenames = {...
    '10987-010_MarkusBauer_20150717';
    'S2-09467-008_MarkusBauer_20150721';
    'S3-11133-01_MarkusBauer_20150724';
    'S4-11139-001_MarkusBauer_20150721';
    'mc09926-004_MarkusBauer_20150722';
    'ia11147-001_MarkusBauer_20150723';
    's8-10862-004_MarkusBauer_20150729';
    'S9-11143-001_MarkusBauer_20150724';
    'S10-11153-001_MarkusBauer_20150731';
    '11151-001_MarkusBauer_20150728';
    'S12-11126-003_MarkusBauer_20150729';
    'S13-11024-007_MarkusBauer_20150804';
    'S14-11013-001_MarkusBauer_20150804';
    's16-11259-001_MarkusBauer_20151009';
    's17-11245-001_MarkusBauer_20151004';
    'S18-11261-001_MarkusBauer_20151009';
    's19-11246-001_MarkusBauer_20151009';
    '11254-001_MarkusBauer_20151007';
    'jr11255-001_MarkusBauer_20151007';
    's22-11270-001_MarkusBauer_20151012';
    's23-11266-001_MarkusBauer_20151012';
    's24-11260-001_MarkusBauer_20151009';
    's25-11258-001_MarkusBauer_20151009';
    's26-11323-001_MarkusBauer_20151030';
    };

sessions = {...
    1;
    [2:5]; %tom
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [2:4];
    [1:3];
    [1:4];
    [2:4];
    [2:4];
    [2:4];
    };

end
