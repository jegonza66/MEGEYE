function fun_rawdatacomponentrejection_v2(datain, datain_eyemap, dataout)
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

    [datdum]        = ft_rejectcomponent(rccfg,comp,datdum);        

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
