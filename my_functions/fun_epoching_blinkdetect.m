function fun_epoching_blinkdetect(icfg)
    cfg.dataset                    = icfg.dataset;
    cfg.continuous                 = 'yes';

    %% define trials / epochs
    cfg.trialdef.eventvalue         = [11:16,21:26,31:36,41:46,51:56,61:66,71:76]+10; %10 times the probab + stim id (face/house) + coherence level
    cfg.trialdef.eventtype          = {'UPPT002'};
    cfg.trialdef.prestim            = 1.3;
    cfg.trialdef.poststim           = 2;
    cfg = ft_definetrial(cfg);
    %cfg.trl([1,end],:) = [];
    save(icfg.preprocfname,'cfg');

    cfg.bsfilter                    = 'yes';
    cfg.bsfreq                      = [49,51];
    cfg.padding                     = 5;
    cfg.blc                         = 'yes';%'no';
    cfg.detrend                     = 'no';
    cfg.artfctdef.feedback          = 'yes';
    cfg.artfctdef.reject            = 'partial';
    cfg.artfctdef.minaccepttim      = 0.1;
    cfg.artfctdef.type              = {'eog'};
    cfg.artfctdef.eog.bpfilter     = 'yes';
    cfg.artfctdef.eog.bpfilttype   = 'but';
    cfg.artfctdef.eog.bpfreq       = [1 15];
    cfg.artfctdef.eog.bpfiltord    = 4;
    cfg.artfctdef.eog.hilbert      = 'yes';
    cfg.artfctdef.eog.channel      = {'EEG057';'EEG058'};
    cfg.artfctdef.eog.cutoff       = 4;   % z-value at which to threshold
    cfg.artfctdef.eog.trlpadding   = 0;   % indicates range to search 4 artifact
    cfg.artfctdef.eog.fltpadding   = 1;   % indicates padding range to avoid filter-artifacts
    cfg.artfctdef.eog.artpadding   = 0.4; % indicates range of data rejected beyond suprathresh.
    cfg.artfctdef.eog.feedback     = 'yes';

    cfgeye = ft_rejectartifact(cfg);
    save(icfg.prepreyefname,'cfgeye');

end