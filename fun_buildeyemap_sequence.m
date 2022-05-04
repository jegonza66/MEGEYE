function sequence = fun_buildeyemap_sequence(cfg)
% function fun_buildeyemap_sequence
% MJI, 01/01/2019
% This function constructs the eyemap sequence following the triggers used in the code. 
% QUESTIONS: Juan, why do we want to change the labels rather than using
% the ones from the code?
% Juan, when adding a condition based on the samples, the samples from
% cfg.trl and the ones here don't seem consistent. Any idea why?

% % 12 Event Marks
%     cfg.evts.BgnEyeMap      = 2; %201;
%     cfg.evts.EndEyeMap      = 4; %202;
%     cfg.evts.OnHL           = 211;
%     cfg.evts.OffHL          = 212;
%     cfg.evts.OnHS           = 221;
%     cfg.evts.OffHS          = 222;
%     cfg.evts.OnVL           = 231;
%     cfg.evts.OffVL          = 232;
%     cfg.evts.OnVS           = 241;
%     cfg.evts.OffVS          = 242;
%     cfg.evts.OnBL           = 251;
%     cfg.evts.OffBL          = 252;

eyemap_triggers = [211,212,221,222,231,232,241,242,251,252];
eyemap_labels   = {'LH','OFF','SH','OFF','LV','OFF','SV','OFF','BL','OFF'}
sequence.equivalence = {'SH', 221; 'LH', 211; 'SV',231; 'LV', 231; 'BL', 251; 'OFF', 200};

for ii=1:size(cfg.trl,1)
    ind=0;
    BgnEyeMap_t(ii)=cfg.trl(ii,1);
    EndEyeMap_t(ii)=cfg.trl(ii,2);
    for ie=1:numel(cfg.event) %loop over all events
        if(isequal(cfg.event(ie).type,string(cfg.trialdef.eventtype))  &&...
                ~isempty(intersect(cfg.event(ie).value,eyemap_triggers)))% &&...
%                 cfg.event(ie).sample>=BgnEyeMap_t(ii) && cfg.event(ie).sample<EndEyeMap_t(ii))
            trigger_value(ii,ind+1)=cfg.event(ie).value;
            trigger_sample(ii,ind+1)=cfg.event(ie).sample;
            [C, ia, ib] = intersect(cfg.event(ie).value,eyemap_triggers);
            sequence.info{ind+1}=eyemap_labels{ib};
            %trigger_value(ii,ind+1)=200;
            %trigger_sample(ii,ind+1)=nan;
            ind=ind+1;
        end
    end
end
sequence.trigger_value=trigger_value;
sequence.trigger_sample=trigger_sample;

%if ~exist([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'file')
%     sequence.info = {'SH', 'OFF', 'SV', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'BL', 'OFF', ...
%         'SH', 'OFF', 'LH', 'OFF', 'BL', 'OFF', 'SV', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF', ...
%         'SV', 'OFF', 'SH', 'OFF', 'BL', 'OFF', 'LH', 'OFF', 'SH', 'OFF', 'SV', 'OFF', 'LV', 'OFF', 'LH', 'OFF', 'LV', 'OFF', 'BL', 'OFF'};
%     sequence.value = [101, 200, 103, 200, 103, 200, 104, 200, 104, 200, 102, 200, 101, 200, 104, 200, 102, 200, 104, 200, ...
%         101, 200, 102, 200, 104, 200, 103, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200, ...
%         103, 200, 101, 200, 104, 200, 102, 200, 101, 200, 103, 200, 104, 200, 102, 200, 104, 200, 104, 200]';
%     sequence.equivalence = {'SH', 101; 'LH', 102; 'SV', 103; 'LV', 104; 'BL', 105; 'OFF', 200};
%    save([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'],'sequence');
%else
%    load([cfg.matdir,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat']);
%end



end