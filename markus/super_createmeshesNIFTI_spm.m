function super_createmeshesNIFTI_spm(mysubjects)

%getdirectories;
experimentdir = 'd:\Projects\MSc_AttentionPrediction\MEG\';
rootmridir = 'd:\Projects\MSc_AttentionPrediction\MR_anatomicals\';
cd(experimentdir)

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
    
    %% try to find a structural MRI in specified directory
    [relentries,entries] = find_directory_entries([rootmridir,directories{isubject}],'nii');
    if ~isempty(relentries)
        %give a consistent name/pattern for available structurals
        %numindx = strmatch([rootmridir,directories{isubject},'\',relentries]); %directories{isubject},'_structural'],
        numindx = 1;
        sMRI = [relentries{numindx},',1'];
    else %if not the default MNI brain will be used subsequently, as indiciated by empty MRI strucutre
        sMRI = [];
    end;
    
    
    %% prepare the forward model - copy mesh into each spm dataset
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        %cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        cd([experimentdir,subjectinfo.directories{(isubject)}]);
        
        spmfile = [cfg.directory,'matfiles\spm8raw_',cfg.myfname([1:(end-3)]),'.mat'];
        
        
        D = spm_eeg_load(spmfile);
        val = 1;
        %[D,val] = spm_eeg_inv_check(D,1,0);
        if ~isfield(D, 'inv') || ~isfield(D.inv{val}, 'comment')
            D.inv = {struct('mesh', [])};
            D.inv{val}.date    = strvcat(date,datestr(now,15));
            D.inv{val}.comment = {''};
        else
            inv = struct('mesh', []);
            inv.comment = D.inv{val}.comment;
            inv.date    = D.inv{val}.date;
            D.inv{val} = inv;
        end
        
        cd([experimentdir,subjectinfo.directories{(isubject)},'\matfiles']);
        
        D.inv{val}.mesh = spm_eeg_inv_mesh(sMRI,2);
        %     D = meeg(D);
        %     %save(spmfile,'D')
        save(D);
        
        relentries = []; entries = [];
        
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
    };

end
