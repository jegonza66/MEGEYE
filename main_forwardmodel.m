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

%% It starts with the conversion of the MEG dataset (‘convert_rawds2spm’)
% super_convert_rawds2SPM(mysubjects)
for isubject = 1:length(directories)
% isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        cfg.dataset = [cfg.directory,cfg.myfname];
        
        if exist(cfg.dataset,'dir')
            spmcfg.dataset          = cfg.dataset; %[sessionfilenames{j},'0',num2str(session_numbers{j}(i)),'.ds',delim,'',sessionfilenames{j},'0',num2str(session_numbers{j}(i)),'.meg4'];
            spmcfg.outfile          = [cfg.directory,'matfiles\spm12raw_',cfg.myfname([1:(end-3)]),'.mat'];
            spmcfg.continuous       = 1;
            spmcfg.checkboundary    = 0;
            spmcfg.channels         = labels;
            spmcfg.inputformat      = [];
            spmcfg.usetrials        = 1;
            spmcfg.timewindow       = [];
            spmcfg.datatype         = 'float32-le'; % (2018-03-27, JK) I get this message: Data type is missing or incorrect, assigning default.
            spmcfg.eventpadding     = 0;
            spmcfg.saveorigheader   = 0;
            spmcfg.conditionlabel   = {'Undefined'};
            D = spm_eeg_convert(spmcfg);
%             clear spmcfg
%             clear D data
        end;
    end;
end;

%% Then the creation of the meshes (if having the anatomical MRI – otherwise it loads the defauilt MNI brain)
% super_createmeshesNIFTI_spm(mysubjects)
for isubject = 1:length(directories)
% isubject = 1;
    
    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    %% try to find a structural MRI in specified directory
%     [relentries,entries] = find_directory_entries([rootmridir,directories{isubject}],'nii');
%     if ~isempty(relentries)
%         %give a consistent name/pattern for available structurals
%         %numindx = strmatch([rootmridir,directories{isubject},'\',relentries]); %directories{isubject},'_structural'],
%         numindx = 1;
%         sMRI = [relentries{numindx},',1'];
%     else %if not the default MNI brain will be used subsequently, as indiciated by empty MRI strucutre
        sMRI = [];
%     end;
        
    %% prepare the forward model - copy mesh into each spm dataset
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        %cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        origpath = pwd;
        cd([experimentdir,subjectinfo.directories{(isubject)}]);
        
        spmfile = [cfg.directory,'matfiles\spm12raw_',cfg.myfname([1:(end-3)]),'.mat'];
        
        
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
        
        cd([experimentdir,subjectinfo.directories{(isubject)},'/matfiles']);
        
        D.inv{val}.mesh = spm_eeg_inv_mesh(sMRI,2);
        %     D = meeg(D);
        %     %save(spmfile,'D')
        save(D);
        
        relentries = []; entries = [];
        
        cd(origpath);
        
    end;
    
end;

%% The coregistration of the mesh surface with the MEG dataset
% super_coregister_spm_v2(mysubjects)
% for isubject = 1:length(directories)
isubject = 1;
    cfg.directory = [experimentdir,directories{isubject},filesep];
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        %cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        origpath = pwd;
        cd([experimentdir,subjectinfo.directories{(isubject)}]);
        
        spmfile = [cfg.directory,'matfiles\spm12raw_',cfg.myfname([1:(end-3)]),'.mat'];
        
        D = spm_eeg_load(spmfile);
        myfiducials = D.inv{1}.mesh.fid;
        myfiducials.fid.pnt = myfiducials.fid.pnt([1,4,5],:);
        %this has to be like that; the labels need to match!
        myfiducials.fid.label = myfiducials.fid.label([1,2,3]); %this has to be like that; the labels need to match!
        
        D = spm_eeg_inv_datareg_ui(D,1,D.fiducials,myfiducials);
        inversion = D.inv;
        save(D);
        clear D;
        clear inversion;

        cd(origpath);
        
    end;
    
% end;

%% The creation of the forward model
% super_createforwardmodel_spm(mysubjects)
% for isubject = 1:length(directories)
isubject = 1;

    cfg.directory = [experimentdir,directories{isubject},filesep];
    
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})
        
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        %cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        origpath = pwd;
        cd([experimentdir,subjectinfo.directories{(isubject)}]);
        
        spmfile = [cfg.directory,'matfiles\spm12raw_',cfg.myfname([1:(end-3)]),'.mat'];
        
        D = spm_eeg_load(spmfile);
        val = 1;
        [D, val] = spm_eeg_inv_check(D,val);
        D.inv{val}.forward = struct([]);
        D.inv{val}.forward(1).voltype = 'Single Shell';
        D = spm_eeg_inv_forward(D);
        spm_eeg_inv_checkforward(D, val);
        
        [L,D] = spm_eeg_lgainmat(D);
        save(D);
        
        fname = [D.path,filesep,'fwdmdl',D.fname];
        forward = D.inv{val};
        save(fname,'forward','L');

        cd(origpath)
    end;
    
% end;

%% Optionally to check that forward model
% super_checkforwardmodels(mysubjects)
% for isubject = 1:length(directories)
isubject = 1;

    cfg.directory = [experimentdir,directories{isubject},filesep];
    display('________________________________________________');
    display(['subject: ',num2str(isubject)]);
    
    for isess = 1:length(sessions{isubject})       
        cfg.myfname = [sessionfilenames{isubject},'_0',num2str(sessions{isubject}(isess)),'.ds'];
        %cfg.dataset = [cfg.directory,filesep,cfg.myfname];
        origpath=pwd;
        cd([experimentdir,subjectinfo.directories{(isubject)}]);
        
        spmfile = [cfg.directory,'matfiles\spm12raw_',cfg.myfname([1:(end-3)]),'.mat'];
        
        D = spm_eeg_load(spmfile);
        val = 1;
        %   Checks that the EEG/MEG .mat file structure is loaded properly and that
        %   the particular inversion of interest has been specified
        [D, val] = spm_eeg_inv_check(D,val); 

        % Check M/EEG forward model
%         The results of coregistration will be presented in SPM’s graphics window. It is important
% to examine the results carefully before proceeding. In the top plot you will see the scalp, the
% inner skull and the cortical mesh with the sensors and the fiducials. For EEG make sure that the
% sensors are on the scalp surface. For MEG check that the head positon in relation to the sensors
% makes sense and the head does not for instance stick outside the sensor array. In the bottom plot
% the sensor labels will be shown in topographical array. Check that the top labels correspond to
% anterior sensors, bottom to posterior, left to left and right to right and also that the labels are
% where you would expect them to be topographically.

        % But, this function is only showing me one plot!!!
        spm_eeg_inv_checkforward(D, val); 
        
%         % Undefined function or variable 'checkforwardmodels_coregistration'.
%         [fighandles] = checkforwardmodels_coregistration(D.inv{1});
%         keyboard;
%         
%         close(fighandles(1));%close(fighandles(2));close(fighandles(3));
        cd(origpath)
    end;
% end;

%%

