function fun_readdata(icfg,useeye)
    if nargin<2
        useeye=0;
    end
    load(icfg.preprocfname,'cfg');
    
    %cfg.channel         = {'MEG'};
    cfg.channel         = {'meg';'megref';'UADC*'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
    cfg.bsfilter        = 'yes';
    cfg.bsfreq          = [49.5,50.5];
    cfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
    cfg.blc             = 'yes';
    cfg.detrend         = 'no';
    data = ft_preprocessing(cfg);

    if (useeye==1)
        save(icfg.datafname,'data');    
        clear data

        load(icfg.prepreyefname,'cfgeye');
        cfg.trl            = cfgeye.artfctdef.eog.artifact;
        cfg.trl(:,3)       = 0;
        mycfg = cfg;
        dateye              = ft_preprocessing(mycfg);
        clear mycfg;
        
        save(icfg.FTeyedatafname,'dateye');
        clear dateye
    else
        ccfg.resamplefs     = 600;      % (JK, 21/03/2018) I leave it like this, but I guess we'll want to match the eye tracking sampling rate
        ccfg.demean         = 'yes';
        ccfg.detrend        = 'yes';    % (JK, 21/03/2018) I leave it like this, but I guess we'll want to change this.
        [data] = ft_resampledata(ccfg, data);

        % (Markus note) Now the data with reference sensors (aligned with the MEG sensors in time)
        % and ocular artefact removed are transformed to, e.g. 1st grad
        gcfg.gradient = 'G3BR'; %'none'
        data = ft_denoise_synthetic(gcfg,data);
        save(icfg.datafname_weyes,'data');

        load('labels_CTF269Nott','labels');
        scfg.channel = labels;
        data = ft_selectdata(scfg,data);
        save(icfg.datafname,'data');    
    end
        
    

end


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