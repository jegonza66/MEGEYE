function fun_replace_eventmarkers(datainfile, dataoutfile, sequencefile)
% Sequence definition for the EyeMap trials of MJI (13229-001)
%   In the future I have to clean the eye data to do it automatically. But
%   for now, I also have problems with the eventmarkers there.
%     sequence.info = {'SH', 'OFF', 'SV', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'BL', 'OFF', ...
%                         'SH', 'OFF', 'LH', 'OFF', 'BL', 'OFF', 'SV', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF', ...
%                         'SV', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF'};
%     sequence.value = [101, 200, 103, 200, 103, 200, 104, 200, 104, 200, 102, 200, 101, 200, 104, 200, 102, 200, 104, 200, ...
%                         101, 200, 102, 200, 104, 200, 103, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200, ...
%                         103, 200, 101, 200, 104, 200, 102, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200]';
%     sequence.equivalence = {'SH', 101; 'LH', 102; 'SV', 103; 'LV', 104; 'BL', 105; 'OFF', 200};
%     save([cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'sequence');


%     dboxpath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/';
%     datapath = [dboxpath 'MEG_recordings/data/'];
%     subjid = '13229_001';
%     eye = load([datapath subjid '/' subjid '.mat']);
%     arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials)
% 
%     sequence.value = arrayfun(@(x) find(strcmp(x.info.trialtype,eye.cfg.tasks)),eye.ExpTrials)';
%     sequence.equivalence = {eye.cfg.tasks{1}, 1; eye.cfg.tasks{2}, 2;eye.cfg.tasks{3}, 3};
%     sequence.info = {};
%     for i=1:length(sequence.value)
%         sequence.info{i} = sequence.equivalence{sequence.value(i),1};
%     end
%     save([cfg.directory,'/matfiles/sequence_', cfg.myfname,'.mat'],'sequence');

%     datainfile  = cfg.eyemap_datafname;
%     dataoutfile = cfg.eyemap_datafname;
%     sequencefile= cfg.eyemap_sequence;
    load(datainfile)
    load(sequencefile)
    % Where's trial information?
    %     data = 
    %          label: {269x1 cell}
    %      trialinfo: [60x1 double]
    %     sampleinfo: [60x2 double]
    %          trial: {1x60 cell}
    %           time: {1x60 cell}
    %            cfg: [1x1 struct]
    
    % data.trialinfo = [2 2 2 ... 2] 
    %   It has the eventmarker value
    % data.sampleinfo  
    %   It has the begin and the end points of the epoch
    
    % The information in cfg only corresponds to the first block, the rest
    % of the information is probably lost in the **merge** step
    %
    % data.cfg.trialdef
    %          eventtype: {'UPPT002'}
    %     eventvalue: 1
    %        prestim: 1.2000
    %       poststim: 4
    % data.cfg.trl
    %       168914      175153       -1440           1
    %       186076      192315       -1440           1
    %       203379      209618       -1440           1
    %       220562      226801       -1440           1
    %       237565      243804       -1440           1
    %       ...         ...         ...             ...
    % data.cfg.event
    %       It has all the events but only for the first block. The rest of
    %       the information is probably lost in the merge step.

    if ( length(data.trialinfo)~=length(sequence.value) ) 
        fprintf('ERROR: Check your sequence definition. The length of the old trial\n');
        fprintf('\t\tinformation is %d, and the new one is %d\n',length(data.trialinfo),length(sequence.value));
        return;
    end
    
    data.trialinfo = sequence.value; %%%[2 2 2]'
    
    
    % cfg.trials      = 1xN, trial indices to keep, can be 'all'. You can use logical indexing, where false(1,N) removes all the trials
    % [data] = ft_selectdata(cfg, data, ...)

    slcfg = [];
    slcfg.trials = (data.trialinfo~=200);
    fprintf('Removing %d trials, and keeping the remaining %d trials...\n',sum(slcfg.trials==0),sum(slcfg.trials==1));
    [data] = ft_selectdata(slcfg, data);

    save(dataoutfile,'data','-v7.3')
    
end