function fun_rawdatacomponentrejection_usingEyeMap(datain, datain_eyemap, dataout)
    datain = cfg.datafname;
    datain_eyemap=cfg.eyemap_datafname_clean;

    load(datain,'data');
    % scfg.channel = {'all', '-EOG'};
    % data = ft_selectdata_new(scfg,data);

    load(datain_eyemap,'rccfg','comp','compno');
    
    rcfg.method         = 'summary';
    rcfg.channel        = 'all';
    rcfg.keepchannel    = 'nan';
    rcfg.feedback       = 'yes';
    datdum              = ft_rejectvisual(rcfg,data);   
    close all;

%     cocfg.blc           = 'yes';
%     cocfg.method        = 'pca';
%     cocfg.channel       = {'all'};
%     cocfg.numcomponent  = length(data.label);
%     comp2               = ft_componentanalysis(cocfg,datdum); 
    
    % I need to project the data into the PCs from the eyemap data
    % and include a field called 'trials'
    [datdum]            = ft_rejectcomponent(rccfg,comp,datdum);        

    rcfg.method         = 'summary';
    rcfg.channel        = 'all';
    rcfg.keepchannel    = 'nan';
    rcfg.feedback       = 'yes';
    datdum              = ft_rejectvisual(rcfg,data);   
    close all;
    
    data    = datdum; 
    clear datdum;
    
    comp    = rmfield(comp,'trial'); 
    compno  = rccfg.component;

    save(dataout,'data','rccfg','comp','compno');
end
