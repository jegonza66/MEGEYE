%% Re-epoch data

function [dataMEG,cfg] = fun_reepoch_data_blinks(cfg,data,eyedata)
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
blicount = nan(length(data.time),2);
clear trial time megevts
for tr = 1:length(data.time);
    blicount(tr) = eyedata(tr).Nblink;
    if ( eyedata(tr).Nblink > 0 )
        for i = 1:eyedata(tr).Nblink
            if ~isnan(eyedata(tr).megevts.blinkOnset(i));
                indt0 = find( data.time{tr} >= eyedata(tr).megevts.blinkOnset(i)/1000,1,'first');
                indt = (indt0+epochLimits_smp(1)):(indt0+epochLimits_smp(2));

                if ( ( min(indt) > 0 ) && ( max(indt) <= length(data.time{tr}) ) )
                    k = k + 1;

                    
                    time{k}  = data.time{tr}(indt) - data.time{tr}(indt0);
                    trial{k} = data.trial{tr}(:,indt);
                    sampleinfo(k,:) = data.sampleinfo(tr,1) + [indt(1) indt(end)];

                    megevts(k).blinkOnset       = eyedata(tr).megevts.blinkOnset(i);
                    megevts(k).blinkDuration      = eyedata(tr).megevts.blinkDuration(i);
                else
                    display('WARNING: There are some blinks outside the trial boundaries')
                end
            end
        end
    end
end
fprintf('Total number of blinks included = %d\n',k)
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
