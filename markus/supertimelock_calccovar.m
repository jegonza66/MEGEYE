
function supertimelock_calccovar(mysubjects)

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
        
        datafname = [cfg.directory,'matfiles',filesep,'data_',num2str(isess),'.mat'];
        load(datafname,'data');
        
        if isfield(data,'label') && (length(data.label) > 274)
            data.label = data.label(1:274);
        end;
        
        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvnofilt_',num2str(isess),'.mat'];
        tmlkcfg = [];
        tmlkcfg.vartrllength = 0;
        tmlkcfg.channel = 'MEG';
        tmlkcfg.lpfilter = 'yes';
        tmlkcfg.lpfreq = 45;
        tmlkcfg.keeptrials = 'no';
        tmlkcfg.covariance = 'yes';
        timelock = ft_timelockanalysis(tmlkcfg,data); 
        save(timelockfile,'timelock');
        
        timelockfile = [cfg.directory,'matfiles',filesep,'tmlklcmvhifrqfilt_',num2str(isess),'.mat'];
        tmlkcfg = [];
        tmlkcfg.vartrllength = 0;
        tmlkcfg.channel = 'MEG';
        tmlkcfg.hpfilter = 'yes';
        tmlkcfg.hpfreq = 40;
        tmlkcfg.keeptrials = 'no';
        tmlkcfg.covariance = 'yes';
        timelock = ft_timelockanalysis(tmlkcfg,data); 
        save(timelockfile,'timelock');
        
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


