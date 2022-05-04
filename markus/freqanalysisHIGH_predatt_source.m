function freqanalysisHIGH_predatt_source(cfg,data)


%% HIGH FREQUENCIES - multitaper
fcfg = [];
fcfg.method = 'mtmconvol';
fcfg.channel = 'all';
fcfg.keeptrials = 'yes';
fcfg.output = 'pow';
fcfg.toi = [-0.8:0.05:1.05];
fcfg.foi = [20:5:130];%[45:0.5:75];
fcfg.t_ftimwin = ones(1,length(fcfg.foi))*0.2;
fcfg.tapsmofrq = ones(1,length(fcfg.foi))*10;
fcfg.pad = 4;
fcfg.keeptrials = 'yes';
fcfg.keeptapers = 'no';
freqhi = ft_freqanalysis(fcfg,data);
clear data;

notdone=1;thresh = 8;niter = 0;
while notdone & (niter < 20)
    [mypowspctrm,outlierindcs] = removeoutliers(freqhi.powspctrm,0,thresh,1);
    %get rows/trials with an unacceptably large artevfact
    dummy = squeeze(sum(sum(sum(mypowspctrm,2),3),4));
    myartfcttrls = (find(isnan(dummy)));
    if (length(myartfcttrls) <= size(freqhi.powspctrm,1)/10)
        notdone = 0;
    else
        thresh = thresh + 2;
    end
    niter = niter + 1;
end;
freqhi.powspctrm = single(mypowspctrm);
clear mypowspctrm;

myfname = [cfg.freqhighfname];
save(myfname,'freqhi');

clear freqhi


end


