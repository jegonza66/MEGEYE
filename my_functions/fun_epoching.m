function fun_epoching(icfg)
    try
        cfg.dataset                    = icfg.dataset;
        cfg.continuous                 = 'yes';

        tmpcfg = cfg;
    %% define trials / epochs
    %     cfg.trialdef.eventtype          = 'gui';

        cfg.trialdef.eventtype          = {'UPPT002'};  % Event Channel
        cfg.trialdef.eventvalue         = [1];          % Values on the event channel
                                                        % name of the triggers in brain vision format
        cfg.trialdef.prestim            = 1.2;   
        cfg.trialdef.poststim           = 4;           %4
        cfg = ft_definetrial(cfg);
        %
        %         % The trial definition "trl" is an Nx3 matrix, N is the number of trials.
        % % The first column contains the sample-indices of the begin of each trial
        % % relative to the begin of the raw data, the second column contains the
        % % sample-indices of the end of each trial, and the third column contains
        % % the offset of the trigger with respect to the trial. An offset of 0
        % % means that the first sample of the trial corresponds to the trigger. A
        % % positive offset indicates that the first sample is later than the trigger,
        % % a negative offset indicates that the trial begins before the trigger.
        % %
        % % The trial definition "trl" can contain additional columns besides the
        % % required three that represend begin, end and offset.
        %  212584      218823       -1440           1
        %  230006      236245       -1440           1
        %  247369      253608       -1440           1
        %  264752      270991       -1440           1
        %cfg.trl([1,end],:) = [];
        save(icfg.preprocfname,'cfg');

        if isfield(icfg,'eyemap_preprocfname')
            cfg = tmpcfg;
            cfg.trialdef.eventtype          = {'UPPT002'};  % Event Channel
            cfg.trialdef.eventvalue         = [2];          % Values on the event channel
                                                            % name of the triggers in brain vision format
            cfg.trialdef.prestim            = 0;   
            cfg.trialdef.poststim           = 6;
%             cfg.em.time     = 6;
%             cfg.em.wait     = 1;
%             cfg.em.sttags   = {'HL','VL','HS','VS','BL','HL','VL','HS','VS','BL'};
            
            cfg = ft_definetrial(cfg);
            %cfg.trl([1,end],:) = [];
            save(icfg.eyemap_preprocfname,'cfg');
        end
    catch ME
        keyboard
    end
end
