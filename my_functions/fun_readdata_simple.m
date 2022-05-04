function fun_readdata_simple(icfg,channels)
    load(icfg.preprocfname,'cfg');
    
    %cfg.channel         = {'MEG'};
    if nargin<2
        cfg.channel         = {'meg';'megref';'UADC*'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
        %cfg.channel         = {'meg';'megref'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
    else 
        cfg.channel         = channels; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
    end
    
    cfg.bsfilter        = 'yes';
    cfg.bsfreq          = [49.5,50.5];
    cfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
    cfg.demean          = 'yes';    % (JK, 22/03/2018) cfg.blc is going to be deprecated
    cfg.detrend         = 'no';
    
    labels_eye={'UADC001';'UADC002';'UADC013'}; %X, Y, Pupil Size

    data = ft_preprocessing(cfg);
    
    ccfg.resamplefs     = 600;      % (JK, 21/03/2018) I leave it like this, but I guess we'll want to match the eye tracking sampling rate
    ccfg.demean         = 'yes';
    ccfg.detrend        = 'no';     % (JK, 21/03/2018) I leave it like this, but I guess we'll want to change this.
    [data_full] = ft_resampledata(ccfg, data);

    ecfg.channel = labels_eye;
    data = ft_selectdata(ecfg,data_full);
    save(icfg.datafname_eyech_only,'data');
    %data_eye = data;
    
    data = data_full;
    % (Markus note) Now the data with reference sensors (aligned with the MEG sensors in time)
    % and ocular artefact removed are transformed to, e.g. 1st grad
    gcfg.gradient = 'G3BR'; %'none'
    data = ft_denoise_synthetic(gcfg,data);
    save(icfg.datafname_weyech,'data');
    %MJI, 28.2.19: By comparing the traces of eye channels before and after
    %gradient correction, we verified that gradients doesn't affect eye
    %channels.
    %figure;plot(data.trial{160}(297,:));hold on;plot(data_eye{160}(1,:),'r--')
    load('labels_CTF269Nott','labels');
    scfg.channel = labels;
    data = ft_selectdata(scfg,data);
    save(icfg.datafname,'data');    
        
    if (isfield(icfg,'eyemap_datafname') && exist(icfg.eyemap_preprocfname))   
        load(icfg.eyemap_preprocfname,'cfg');
    
        %cfg.channel         = {'MEG'};
        cfg.channel         = {'meg';'megref';'UADC*'}; % (JK, 21/03/2018) I include the two eye tracking channels (UADC013 and UADC014)
        cfg.bsfilter        = 'yes';
        cfg.bsfreq          = [49.5,50.5];
        cfg.padding         = 7.168;    % (JK, 22/03/2018) How is calculated? 
        cfg.demean          = 'yes';    % (JK, 22/03/2018) cfg.blc is going to be deprecated
        cfg.detrend         = 'no';
        data = ft_preprocessing(cfg);

        ccfg.resamplefs     = 600;      % (JK, 21/03/2018) I leave it like this, but I guess we'll want to match the eye tracking sampling rate
        ccfg.demean         = 'yes';
        ccfg.detrend        = 'yes';    % (JK, 21/03/2018) I leave it like this, but I guess we'll want to change this.
        [data_full] = ft_resampledata(ccfg, data);

        ecfg.channel = labels_eye;
        data = ft_selectdata(ecfg,data_full);
        save(icfg.eyemap_datafname_eyech_only,'data');
        data = data_full;
        % (Markus note) Now the data with reference sensors (aligned with the MEG sensors in time)
        % and ocular artefact removed are transformed to, e.g. 1st grad
        gcfg.gradient = 'G3BR'; %'none'
        data = ft_denoise_synthetic(gcfg,data);
        save(icfg.eyemap_datafname_weyech,'data');

        load('labels_CTF269Nott','labels');
        scfg.channel = labels;
        data = ft_selectdata(scfg,data);
        save(icfg.eyemap_datafname,'data');    
    end    

end