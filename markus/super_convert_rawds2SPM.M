
function super_convert_rawds2SPM(mysubjects)

experimentdir = 'd:\Projects\MSc_AttentionPrediction\MEG\';
[directories,sessionfilenames,sessions] = collectsubjectinfo;
load(['labels_CTF275.mat']);

if nargin<1
    mysubjects = [1:length(directories)];
end;
mysubjects


for isubject = mysubjects
    
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        cfg.dataset = [cfg.directory,cfg.myfname];
        
        
        if exist(cfg.dataset,'dir')
            
            spmcfg.dataset = cfg.dataset; %[sessionfilenames{j},'0',num2str(session_numbers{j}(i)),'.ds',delim,'',sessionfilenames{j},'0',num2str(session_numbers{j}(i)),'.meg4'];
            spmcfg.outfile = [cfg.directory,'matfiles\spm8raw_',cfg.myfname([1:(end-3)]),'.mat'];
            spmcfg.continuous = 1;
            spmcfg.checkboundary = 0;
            spmcfg.channels = labels;
            spmcfg.inputformat = [];
            spmcfg.usetrials = 1;
            spmcfg.timewindow = [5,100];
            spmcfg.datatype = 'float32-le';
            spmcfg.eventpadding = 0;
            spmcfg.saveorigheader = 0;
            spmcfg.conditionlabel = {'Undefined'};
            D = spm_eeg_convert(spmcfg);
            clear spmcfg
            clear D data
            
        end;
        
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
