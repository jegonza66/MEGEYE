function fun_rawdatacomponentrejection(datain, dataout)
    load(datain,'data');

    % scfg.channel = {'all', '-EOG'};
    % data = ft_selectdata_new(scfg,data);

    startover = 1;
    while startover
        rcfg.method         = 'summary';
        rcfg.channel        = 'all';
        rcfg.keepchannel    = 'nan';
        rcfg.feedback       = 'yes';
        datdum              = ft_rejectvisual(rcfg,data);   

        cocfg.blc           = 'yes';
        cocfg.method        = 'pca';
        cocfg.channel       = {'all'};
        cocfg.numcomponent  = length(data.label);
        comp                = ft_componentanalysis(cocfg,datdum);
        % unmixing*data
        
        display('PCA of raw temp data');
        pltcfg.layout = 'CTF269Nott.lay';
        pltcfg.marker = 'off';
        pltcfg.comment= 'no';
        figure; set(gcf,'Position',[675 50 1225 900])
            for pp = 1:8
                pltcfg.component = pp;
                subplot(2,4,pp);
                    if size(comp.topo,1) == 269
                        ft_topoplotIC(pltcfg, comp)
%                         topoplot(comp.topo(:,pp),'layout','CTF269Nott.lay','electrodes','off');
                    elseif size(comp.topo,1) == 270
                        ft_topoplotIC(pltcfg, comp)
%                         topoplot(comp.topo(:,pp),'layout','CTF270Nott.lay','electrodes','off');
                    end;
            end;   
        
        rccfg.component = input('Which bad PCA components do you want to remove? E.g., 1 or [1,2] or [] ');
        [datdum]        = ft_rejectcomponent(rccfg,comp,datdum);        

        rcfg.method     = 'summary';
        rcfg.channel    = 'all';
        rcfg.keepchannel= 'nan';
        rcfg.feedback   = 'yes';
        datdum          = ft_rejectvisual(rcfg,datdum);   

        startover = input('Do you want to start all over again with trial rejection ? [0=no / 1=yes] ');
        clear cocfg;
        close all;
    end;

    data    = datdum; 
    clear datdum;
    comp    = rmfield(comp,'trial'); 
    compno  = rccfg.component;

    save(dataout,'data','rccfg','comp','compno');

end
