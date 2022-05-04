function super_freqanalysis_source_ROISAM_v2hires(mysubjects)

experimentdir = 'd:\Projects\MSc_AttentionPrediction\MEG\';
cd(experimentdir)

mysurfhi = load([experimentdir,'analysisroutines\canonical_surf_MNI.mat']);
mysurflo = load([experimentdir,'analysisroutines\canonical_surf_MNI_413.mat']);


[directories,sessionfilenames,sessions] = collectsubjectinfo;
subjectinfo.directories = directories;
subjectinfo.sessionfilenames = sessionfilenames;

if nargin<1
    mysubjects = [1:length(directories)];
end;
mysubjects

for isubject = mysubjects
    
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        datafname = [cfg.directory,'matfiles',filesep,'data_',num2str(isess),'.mat'];
        load(datafname,'data');
        
        data.label = data.label(1:274);
        if (size(data.trial{1},1) > 274) & (length(data.label) <= 274)
            for itrl=1:length(data.trial)
                data.trial{itrl} = data.trial{itrl}([1:length(data.label)],:);
            end;
        end;
        
        scfg.latency = [-1.2,1.19];
        scfg.channel = {'MEG'}; %;'-MLO12';'-MLO23';'-MRO41';'-MRO42'};
        data = ft_selectdata(scfg,data);
        
        %removal of flicker data
        data = projector_flicker_removal_bandstop(data);
                
        %% LOW frequencies
        filterfile = [cfg.directory,'matfiles',filesep,'LCMVBNDEMnofiltSAMhires_',num2str(sessions{isubject}(isess))];
        load(filterfile);
        
        searchregion = 7; %7; %mm
        sconfig = struct();
        sconfig.filter = findbeamfcoeffs4ROIlowreshiresgrid(mysurflo.surf,mysurfhi.surf,source,searchregion);
        sconfig.sourcelabel = source.cfg.previous.channel;
        sconfig.preproc.lpfilter = 'yes'; 
        sconfig.preproc.lpfreq = 40;        
        clear filter source
        [datdum] = projectdatathrough_allfilters_ROIPCAlatest(data,sconfig);
        
        cfg.freqlowfname = [cfg.directory,'matfiles',filesep,'freqLowhiressess',num2str(isess)];
        freqanalysisLOW_predatt_source(cfg,datdum);
        clear datdum;
        
        
        %% HIGH frequencies
        filterfile = [cfg.directory,'matfiles',filesep,'LCMVBNDEMhifrqfiltSAMhires_',num2str(sessions{isubject}(isess))];
        load(filterfile);
        
        searchregion = 7; %mm
        sconfig = struct();
        sconfig.filter = findbeamfcoeffs4ROIlowreshiresgrid(mysurflo.surf,mysurfhi.surf,source,searchregion);
        sconfig.sourcelabel = source.cfg.previous.channel;
        sconfig.preproc.hpfilter = 'yes'; 
        sconfig.preproc.hpfreq = 35;
        clear filter source
        [datdum] = projectdatathrough_allfilters_ROIPCAlatest(data,sconfig);
        
        cfg.freqhighfname = [cfg.directory,'matfiles',filesep,'freqHighhiressess',num2str(isess)];
        freqanalysisHIGH_predatt_source(cfg,datdum);
        clear source datdum data;
        
    end;
    
end;

end


function freqanalysisLOW_predatt_source(cfg,data)

%% LOW FREQUENCIES
fcfg = [];
fcfg.method = 'mtmconvol';
fcfg.channel = 'all';
fcfg.keeptrials = 'yes';
fcfg.taper = 'hanning';
fcfg.output = 'pow';
fcfg.toi = [-1:0.1:0.8];
fcfg.foi = [2:1:40];%[45:0.5:75];
fcfg.t_ftimwin = ones(1,length(fcfg.foi))*0.5;
fcfg.pad = 4;
fcfg.keeptrials = 'yes';
freqlow = ft_freqanalysis(fcfg,data);

% dummy = (freqlow.powspctrm(:,1:413,:,:) + freqlow.powspctrm(:,414:826,:,:))/2;
% freqlow.powspctrm = dummy;
% freqlow.label = freqlow.label(1:413);


notdone=1;thresh = 8;niter = 0;
while notdone & (niter < 20)
    [mypowspctrm,outlierindcs] = removeoutliers(freqlow.powspctrm,0,thresh,1);
    %get rows/trials with an unacceptably large artevfact
    dummy = squeeze(sum(sum(sum(mypowspctrm,2),3),4));
    myartfcttrls = (find(isnan(dummy)));
    if (length(myartfcttrls) <= size(freqlow.powspctrm,1)/10)
        notdone = 0;
    else
        thresh = thresh + 2;
    end
    niter = niter + 1;
end;
freqlow.powspctrm = single(mypowspctrm);
%update the trialinfo!!!!!
clear mypowspctrm;

myfname = [cfg.freqlowfname];
save(myfname,'freqlow');

clear freqlow;

end


function freqanalysisHIGH_predatt_source(cfg,data)


%% HIGH FREQUENCIES - multitaper
fcfg = [];
fcfg.method = 'mtmconvol';
fcfg.channel = 'all';
fcfg.keeptrials = 'yes';
fcfg.output = 'pow';
fcfg.toi = [-0.8:0.05:1.05];
fcfg.foi = [20:5:130];%[45:0.5:75];
fcfg.t_ftimwin = ones(1,length(fcfg.foi))*0.2;
fcfg.tapsmofrq = ones(1,length(fcfg.foi))*10;
fcfg.pad = 4;
fcfg.keeptrials = 'yes';
fcfg.keeptapers = 'no';
freqhi = ft_freqanalysis(fcfg,data);
clear data;

notdone=1;thresh = 8;niter = 0;
while notdone & (niter < 20)
    [mypowspctrm,outlierindcs] = removeoutliers(freqhi.powspctrm,0,thresh,1);
    %get rows/trials with an unacceptably large artevfact
    dummy = squeeze(sum(sum(sum(mypowspctrm,2),3),4));
    myartfcttrls = (find(isnan(dummy)));
    if (length(myartfcttrls) <= size(freqhi.powspctrm,1)/10)
        notdone = 0;
    else
        thresh = thresh + 2;
    end
    niter = niter + 1;
end;
freqhi.powspctrm = single(mypowspctrm);
clear mypowspctrm;

myfname = [cfg.freqhighfname];
save(myfname,'freqhi');

clear freqhi


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


