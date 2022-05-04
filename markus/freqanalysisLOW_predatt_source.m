function freqanalysisLOW_predatt_source(cfg,data)

%% LOW FREQUENCIES
fcfg = [];
fcfg.method         = 'mtmconvol';
fcfg.channel        = 'all';
fcfg.keeptrials     = 'yes';
fcfg.taper          = 'hanning';
fcfg.output         = 'pow';
fcfg.toi            = [-1:0.1:0.8];
fcfg.foi            = [2:1:40];%[45:0.5:75];
fcfg.t_ftimwin      = ones(1,length(fcfg.foi))*0.5;
fcfg.pad            = 4;
fcfg.keeptrials     = 'yes';
freqlow = ft_freqanalysis(fcfg,data);

% dummy = (freqlow.powspctrm(:,1:413,:,:) + freqlow.powspctrm(:,414:826,:,:))/2;
% freqlow.powspctrm = dummy;
% freqlow.label = freqlow.label(1:413);


notdone=1;thresh = 8;niter = 0;
while notdone & (niter < 20)
    [mypowspctrm,outlierindcs] = removeoutliers(freqlow.powspctrm,0,thresh,1);
    %get rows/trials with an unacceptably large artevfact
    dummy = squeeze(sum(sum(sum(mypowspctrm,2),3),4));
    myartfcttrls = (find(isnan(dummy)));
    if (length(myartfcttrls) <= size(freqlow.powspctrm,1)/10)
        notdone = 0;
    else
        thresh = thresh + 2;
    end
    niter = niter + 1;
end;
freqlow.powspctrm = single(mypowspctrm);
%update the trialinfo!!!!!
clear mypowspctrm;

myfname = [cfg.freqlowfname];
save(myfname,'freqlow');

clear freqlow;

end


