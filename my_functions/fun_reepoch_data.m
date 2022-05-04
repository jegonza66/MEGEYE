%% Re-epoch data

function [dataMEG,cfg] = fun_reepoch_data(cfg,data,eyedata)
% data = 
%          label: {374x1 cell}
%      trialinfo: [160x1 double]
%     sampleinfo: [160x2 double]
%          trial: {1x160 cell}
%           time: {1x160 cell}
%            cfg: [1x1 struct]
%        fsample: 600
% In data doesn't seem to be an event variable to change and re-epoch the
% data. But, on the other hand, I can do it without changing much.

if ~isfield(cfg,'epochLimits_ms')
    cfg.epochLimits_ms  = [-200 250]; % In ms,
end
epochLimits_smp = round(data.fsample*cfg.epochLimits_ms/1000); % In ms,

if ~isfield(cfg,'thrDuration')
    cfg.thrDuration = [cfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
end

%%
k = 0;
sampleinfo = nan(1000,2);
fixcount = nan(length(data.time),2);
clear trial time megevts
for tr = 1:length(data.time);
    fixcount(tr,1) = eyedata(tr).Nfix;
    if ( eyedata(tr).Nfix > 0 )
        if (cfg.filterNonStimFix)
            filtro = eyedata(tr).megevts.fixType~=0;
        else
            filtro = ones(size(eyedata(tr).megevts.fixType));
        end
        if ~any(isnan(cfg.thrDuration))
            ind = find(filtro & eyedata(tr).megevts.fixDuration >= cfg.thrDuration(1) & eyedata(tr).megevts.fixDuration <= cfg.thrDuration(2));
        elseif isnan(cfg.thrDuration(1))
            ind = find(filtro & eyedata(tr).megevts.fixDuration <= cfg.thrDuration(2));
        elseif isnan(cfg.thrDuration(2))
            ind = find(filtro & eyedata(tr).megevts.fixDuration >= cfg.thrDuration(1));
        else
            ind = find(filtro);
        end

        fixcount(tr,2) = sum(~isnan(eyedata(tr).megevts.fixOnset));
        for i = 1:length(ind)
            if ~isnan(eyedata(tr).megevts.fixOnset(ind(i)));
                indt0 = find( data.time{tr} >= eyedata(tr).megevts.fixOnset(ind(i))/1000,1,'first');
                indt = (indt0+epochLimits_smp(1)):(indt0+epochLimits_smp(2));
                
                if ( ( min(indt) > 0 ) & ( max(indt) <= length(data.time{tr}) ) )
                    k = k + 1;

                    time{k}  = data.time{tr}(indt) - data.time{tr}(indt0);
                    trial{k} = data.trial{tr}(:,indt);
                    sampleinfo(k,:) = data.sampleinfo(tr,1) + [indt(1) indt(end)];

                    megevts(k).fixOnset         = eyedata(tr).megevts.fixOnset(ind(i));
                    megevts(k).fixDuration      = eyedata(tr).megevts.fixDuration(ind(i));
                    megevts(k).fixPosition      = eyedata(tr).megevts.fixPosition(ind(i));
                    megevts(k).fixRank          = eyedata(tr).megevts.fixRank(ind(i));
                    megevts(k).prevfixDuration  = eyedata(tr).megevts.prevfixDuration(ind(i));
                    megevts(k).prevSaccDuration = eyedata(tr).megevts.prevSaccDuration(ind(i));
                    megevts(k).prevSaccAmplitude= eyedata(tr).megevts.prevSaccAmplitude(ind(i));
                    megevts(k).prevSaccDirection= eyedata(tr).megevts.prevSaccDirection(ind(i));
                    megevts(k).fixTrialType     = eyedata(tr).megevts.fixTrialType(ind(i));
                    megevts(k).fixStimType      = eyedata(tr).megevts.fixStimType(ind(i));
                    megevts(k).fixType          = eyedata(tr).megevts.fixType(ind(i));
                    megevts(k).Nfix             = eyedata(tr).Nfix;
                else
                    display('WARNING: There are some fixations outside the trial boundaries')
                end
            end
        end
    end
end
fprintf('Total number of fixations included = %d\n',k)
sampleinfo = sampleinfo(sum(isnan(sampleinfo)')'==0,:);

cfg.megevts = megevts;

%%
% data = 
%          label: {374x1 cell}
%      trialinfo: [160x1 double]
%     sampleinfo: [160x2 double]
%          trial: {1x160 cell}
%           time: {1x160 cell}
%            cfg: [1x1 struct]
%        fsample: 600

    dataMEG = data;
    dataMEG.trialinfo  = ones(size(time))'; % My variable of interest
    dataMEG.sampleinfo = sampleinfo;
    dataMEG.trial      = trial;
    dataMEG.time       = time;
    dataMEG.cfg.trials = ones(size(time))'; % My variable of interest
end
