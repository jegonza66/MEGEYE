load tmpcomp.mat

load(cfg.eyemap_datafname_weyech)
eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/'; 
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata_EyeMap.mat']);
eyedata = eyedata_EyeMap;

eyedata = fun_estimating_lag(data,eyedata,0);
fprintf('Total fixations = %d\n',sum([eyedata.Nfix]))
fprintf('Total blinks = %d\n',sum([eyedata.Nblink]))

% [comp,compno]=fun_rawdatacomponentrejection_v3(cfg.eyemap_datafname,'tmp','ica',{100});

recfg = [];
recfg.filterNonStimFix = 0;
recfg.epochLimits_ms= [-500 500]; % In ms,
recfg.thrDuration   = [50 NaN]; % In ms, [200 NaN] means only larger than 200ms
% recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[compEpoched,recfg] = fun_reepoch_data(recfg,comp,eyedata);

%%
t = compEpoched.time{1};
prevSaccDuration    = [recfg.megevts.prevSaccDuration]/1000;
fixDuration         = [recfg.megevts.fixDuration]/1000;
var_fix     = nan(length(compEpoched.label),length(compEpoched.trial));
var_sacc    = nan(length(compEpoched.label),length(compEpoched.trial));
for ic=1:length(compEpoched.label)
    for tr=1:length(compEpoched.trial)
       var_fix(ic,tr)   = var(compEpoched.trial{tr}( ic , t>0 & t<min([fixDuration(tr),0.5]) ) );
       var_sacc(ic,tr)  = var(compEpoched.trial{tr}( ic , t>(-prevSaccDuration(tr)-0.005) & t<0.010 ) );
    end
end
ratio   = mean(var_sacc./var_fix,2);

thr     = 0.1:0.05:2;
cratio  = nan(size(thr));
for i=1:length(thr); 
    cratio(i) = sum(ratio>thr(i))/length(ratio);
end
figure(1); clf; set(gcf,'Color','w') 
    hold on
        plot([1.1 1.1],[0 1],'k--')
        plot(thr,cratio,'k-','LineWidth',3)
    hold off
    box on
    xlabel('Threshold')
    ylabel('P(var_{sacc}/var_{fix} > Threshold)')

%%
pltcfg.layout   = 'CTF269Nott.lay';
pltcfg.marker   = 'off';
pltcfg.comment  = 'no';
pltcfg.title    = 'off';

[~,ind] = sortrows(ratio,-1);
theta       = [recfg.megevts.prevSaccDirection];
horiz_r     = ( (theta>(2*pi-pi/6)) | (theta<pi/6) );
horiz_l     = ( theta>(pi-pi/6) & theta<(pi+pi/6) );
vert_up     = ( theta>(pi/2 - pi/6) & theta<(pi/2 + pi/6) );
vert_do     = ( theta>(3*pi/2-pi/6) & theta<(3*pi/2+pi/6) );
YLIMI = [-2e-12 2e-12];
icmax = 8;
figure(2); clf; set(gcf,'Color','w','Position',[66 1 1855 1001])
    for ic = 1:icmax
        subplot(icmax,5,5*(ic-1) + 1)
            pltcfg.component = ind(ic);
            ft_topoplotIC(pltcfg, comp)
        
        subplot(icmax,5,5*(ic-1) + 2)
            m=fun_meancurve_cells(compEpoched.trial(vert_up),ind(ic));
            mSaccDur = mean(prevSaccDuration(vert_up));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Upward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); ylabel('Amplitude (fT)'); end
        subplot(icmax,5,5*(ic-1) + 3)
            m=fun_meancurve_cells(compEpoched.trial(vert_do),ind(ic));
            mSaccDur = mean(prevSaccDuration(vert_do));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Downward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
        subplot(icmax,5,5*(ic-1) + 4)
            m=fun_meancurve_cells(compEpoched.trial(horiz_l),ind(ic));
            mSaccDur = mean(prevSaccDuration(horiz_l));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Leftward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
        subplot(icmax,5,5*(ic-1) + 5)
            m=fun_meancurve_cells(compEpoched.trial(horiz_r),ind(ic));
            mSaccDur = mean(prevSaccDuration(horiz_r));
            hold on
                plot(-[mSaccDur mSaccDur],YLIMI,'g--')
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Rightward'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); end
    end
    
%%
blicfg = [];
blicfg.filterNonStim = 0;
blicfg.epochLimits_ms= [-1000 1000]; % In ms,
blicfg.thrDuration   = [50 NaN]; % In ms, [200 NaN] means only larger than 200ms
% recfg.thrDuration   = [recfg.epochLimits_ms(2) 1000]; % In ms, [200 NaN] means only larger than 200ms
[compEpochedbli,blicfg] = fun_reepoch_data_blinks(blicfg,comp,eyedata);

pltcfg.layout   = 'CTF269Nott.lay';
pltcfg.marker   = 'off';
pltcfg.comment  = 'no';
pltcfg.title    = 'off';

YLIMI = [-2e-12 2e-12];
icmax = 8;
figure(3); clf; set(gcf,'Color','w','Position',[66 1 1855 1001])
    for ic = 1:icmax*3
        subplot(icmax,6,2*(ic-1) + 1)
            pltcfg.component = ic;
            ft_topoplotIC(pltcfg, comp)
        
        subplot(icmax,6,2*(ic-1) + 2)
            m=fun_meancurve_cells(compEpoched.trial,ic);
            hold on
                plot([0 0],YLIMI,'r--')
                plot(t,m,'k-')
            hold off
            ylim(YLIMI)
            if (ic==1); title('Blink'); end
            if (ic<icmax); set(gca,'XColor','w'); else xlabel('Time (s)'); ylabel('Amplitude (fT)'); end
    end
