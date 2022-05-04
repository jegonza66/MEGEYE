function [comp,compno] = fun_rawdatacomponentrejection_v3(datain, dataout, method, params)
load(datain,'data');

    % scfg.channel = {'all', '-EOG'};
    % data = ft_selectdata_new(scfg,data);
    switch method
        case 'pca'
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

        case 'ica'
            if (nargin < 4); 
                pcs = input('How many PCA components do you want to keep? E.g., 1 or [] '); 
            else
                pcs = params{1};
            end
            rcfg.method         = 'summary';
            rcfg.channel        = 'all';
            rcfg.keepchannel    = 'nan';
            rcfg.feedback       = 'yes';
            datdum              = ft_rejectvisual(rcfg,data);   

    %         cocfg.blc           = 'yes';
    %         cocfg.method        = 'pca';
    %         cocfg.channel       = {'all'};
    %         cocfg.numcomponent  = length(data.label);
    %         comp                = ft_componentanalysis(cocfg,datdum);
    % 
    %         rccfg.component = (pcs+1):length(comp.label);
    %         [datdum]        = ft_rejectcomponent(rccfg,comp,datdum);        
    %         
            cocfg.binica.extended = 1;
            cocfg.binica.pca  = pcs;
    %         cocfg.binica.sphering
    %         cocfg.binica.lrate
    %         cocfg.binica.blocksize
    %         cocfg.binica.maxsteps
    %         cocfg.binica.stop
    %         cocfg.binica.weightsin
    %         cocfg.binica.verbose
    %         cocfg.binica.filenum
    %         cocfg.binica.posact
    %         cocfg.binica.annealstep
    %         cocfg.binica.annealdeg
    %         cocfg.binica.bias
    %         cocfg.binica.momentum

            cocfg.blc           = 'yes';
            cocfg.method        = 'binica';
            cocfg.channel       = {'all'};
            comp                = ft_componentanalysis(cocfg,datdum);

            fprintf('WARNING: Now Im going to applied the criteria outside the function.\n')
            fprintf('WARNING: Im not substracting the components, nor savig the resultas, only a tmp file.\n')
            save('tmpcomp','comp')
            
%             display('ICA of raw temp data');
%             pltcfg.layout = 'CTF269Nott.lay';
%             pltcfg.marker = 'off';
%             pltcfg.comment= 'no';
%             for i=1:5
%                 figure; set(gcf,'Position',[675 50 1225 900])
%                     for pp = 1:20
%                         pltcfg.component = 20*(i-1)+pp;
%                         subplot(4,5,pp);
%                             if size(comp.topo,1) == 269
%                                 ft_topoplotIC(pltcfg, comp)
%         %                         topoplot(comp.topo(:,pp),'layout','CTF269Nott.lay','electrodes','off');
%                             elseif size(comp.topo,1) == 270
%                                 ft_topoplotIC(pltcfg, comp)
%         %                         topoplot(comp.topo(:,pp),'layout','CTF270Nott.lay','electrodes','off');
%                             end;
%                     end;   
%             end
% 
%             rccfg.component = input('Which bad ICA components do you want to remove? E.g., 1 or [1,2] or [] ');
%             [datdum]        = ft_rejectcomponent(rccfg,comp,datdum);        
% 
%             rcfg.method     = 'summary';
%             rcfg.channel    = 'all';
%             rcfg.keepchannel= 'nan';
%             rcfg.feedback   = 'yes';
%             datdum          = ft_rejectvisual(rcfg,datdum);   
% 
%             data    = datdum; 
%             clear datdum;
%             comp    = rmfield(comp,'trial'); 
%             compno  = rccfg.component;
% 
%             save(dataout,'data','rccfg','comp','compno');

        case 'previous_unmixing'
            comp = params{1};
            compno = params{2};

            rcfg.method         = 'summary';
            rcfg.channel        = 'all';
            rcfg.keepchannel    = 'nan';
            rcfg.feedback       = 'yes';
            datdum              = ft_rejectvisual(rcfg,data);   

    %   Instead of specifying a component analysis method, you can also specify
    %   a previously computed unmixing matrix, which will be used to estimate the
    %   component timecourses in this data. This requires        
            cocfg.unmixing      = comp.unmixing;
            cocfg.topolabel     = comp.topolabel;
            comp                = ft_componentanalysis(cocfg,datdum);

            rccfg.component = compno;
            [datdum]        = ft_rejectcomponent(rccfg,comp,datdum);        

            clear cocfg;
            close all;

            rcfg.method     = 'summary';
            rcfg.channel    = 'all';
            rcfg.keepchannel= 'nan';
            rcfg.feedback   = 'yes';
            datdum          = ft_rejectvisual(rcfg,datdum);   

            data    = datdum; 
            clear datdum;
            comp    = rmfield(comp,'trial'); 
            compno  = rccfg.component;

            save(dataout,'data','rccfg','comp','compno');

    end
end
    
