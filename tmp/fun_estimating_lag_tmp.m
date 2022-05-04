%% 1. How to synchronize the MEG time with the ET time?
clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath('/home/juank/toolbox/fieldtrip-20170827')
% ft_defaults;
% addpath('my_functions/')
experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319_merged.ds'};
                
j = 1;
    cfg.directory   = [experimentdir,directories{j}];
    cfg.myfname     = sessionfilenames{j}([1:(end-3)]);
    cfg.dataset     = [cfg.directory,filesep,sessionfilenames{j}];
    clear experimentdir directories sessionfilenames
    
    display('________________________________________________');
    display(['subject: ',num2str(j) ', ' cfg.myfname]);
    cfg.datafname               = [cfg.directory,'/matfiles/data_', cfg.myfname,'.mat'];
    cfg.datafname_weyech        = [cfg.directory,'/matfiles/data_weyech_', cfg.myfname,'.mat'];
    cfg.preprocfname            = [cfg.directory,'/matfiles/cfg_',  cfg.myfname,'.mat'];
        
    cfg.eyemap_datafname        = [cfg.directory,'/matfiles/eyemap_data_', cfg.myfname,'.mat'];
    cfg.eyemap_datafname_weyech = [cfg.directory,'/matfiles/eyemap_data_weyech_', cfg.myfname,'.mat'];
    cfg.eyemap_preprocfname     = [cfg.directory,'/matfiles/eyemap_cfg_',  cfg.myfname,'.mat'];
     
    cfg.sequence                = [cfg.directory,'/matfiles/sequence_', cfg.myfname,'.mat'];
    cfg.eyemap_sequence         = [cfg.directory,'/matfiles/eyemap_sequence_', cfg.myfname,'.mat'];

%     load(cfg.eyemap_datafname_weyech)
% 
%     indch = find(cellfun(@(x) ~isempty(strfind(x,'UADC')),data.label));
%     data.label(indch);
% 
%     eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
%     subjid = '13229_001';
%     load([eyepath subjid '/' subjid '_eyedata_EyeMap.mat']);
%     eyedata = eyedata_EyeMap;

    load(cfg.datafname_weyech)

    indch = find(cellfun(@(x) ~isempty(strfind(x,'UADC')),data.label));
    data.label(indch);

    eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
    subjid = '13229_001';
    load([eyepath subjid '/' subjid '_eyedata.mat']);
    
%%
xc_xy       = nan(length(data.time),2);
label_xy    = {'Horizontal Position (pxls)','Vertical Position (pxls)'};
for tr = 1:length(data.time);
%     figure(tr); clf
    figure(1); clf
        set(gcf,'Color','w')
        set(gcf,'Position',[675 135 1085 835])
        for i=1:2
            % Analog Channel MEG data
%             t1 = data.time{tr};
%             z1 = data.trial{tr}(indch(i),:);
            t1 = data.time{tr}(data.time{tr}>=0);
            z1 = data.trial{tr}(indch(i),data.time{tr}>=0);
            % Eye Tracking data
            t2 = eyedata(tr).samples(:,1)/1000;
            z2 = eyedata(tr).samples(:,i+1);
            % Interpolate the Eye Tracking data into the MEG timepoints
            z2i = interp1(t2,z2,t1);
            
            % I'm going to interpolate the NaN values with a rect only for
            % the lag estimation
            if any(isnan(z2i))
                indini = find(diff(isnan(z2i))==1);
                indend = find(diff(isnan(z2i))==-1);
                if ( isempty(indini) );  indini = 1;  end
                if ( isempty(indend) );  indend = length(z2i)-1;  end

                if ( indini(1)  > indend(1) );  indini = [1 indini];            end
                if ( indini(end)> indend(end) );indend = [indend length(z2i)-1];end

                for j = 1:length(indini)
                    ti = t1(indini(j));
                    tf = z2i(indini(j));
                    yi = t1(indend(j)+1);
                    yf = z2i(indend(j)+1);
                    
                    if isnan(yi); yi=yf; end
                    if isnan(yf); yf=yi; end

                    a = (yf-yi)/(tf-ti);
                    b = yf - a*tf;

                    z2i(indini(j):(indend(j)+1)) = a*t1(indini(j):(indend(j)+1)) + b;
                end
            end
            
            % Calculate the maximum of the cross-correlation 
            [xcf,lags] = xcorr( (z2i-mean(z2i))/std(z2i) , (z1-mean(z1))/std(z1));
            tlags = 1000*lags/data.fsample;
            [mxcf,imxcf]=max(xcf);
            xc_xy(tr,i) = tlags(imxcf);
            
            subplot(2,2,2*(i-1)+1)
                hold on
                    plot(t1,z1,'b-')
                hold off
                set(gca,'XLim',[min([t1(1) t2(1)]) max([t1(end) t2(end)])])
                ax1 = gca; % current axes
                ax1.XColor = 'b';
                ax1.YColor = 'b';
                ylabel(ax1,'Amplitude')
                xlabel(ax1,'MEG Time (s)')

            ax1_pos = ax1.Position; % position of first axes
            ax2 = axes('Position',ax1_pos,...
                'XAxisLocation','top',...
                'YAxisLocation','right',...
                'Color','none');
                hold on
%                     plot(t2,z2,'r-')
                    plot(t1,z2i,'r-')
                hold off
                set(gca,'XLim',[min([t1(1) t2(1)]) max([t1(end) t2(end)])])
                ax2.XColor = 'r';
                ax2.YColor = 'r';
                ylabel(ax2,'Horizontal position (pxls)')
                xlabel(ax2,'ET Time (s)')
            
            ylabel(ax2,label_xy{i})
            xlabel(ax2,'ET Time (s)')

        subplot(2,2,2*(i-1)+2)
            hold on
                plot(tlags,xcf,'k.-');
                plot(tlags(imxcf),xcf(imxcf),'ro');
                title(tlags(imxcf))
            hold off
            xlabel('lags [ms]'), ylabel('corr')
            xlim([-100 100])
        end
    pause(0.1)
end

%% 
% * Correct samples by a fixed time.
% * Calculate fixation/saccades from samples.
% * Change fixation times to MEG time points.

% xc_xy       = nan(length(data.time),2);
label_xy    = {'Horizontal Position (pxls)','Vertical Position (pxls)'};
for tr = 1:length(data.time);
%     figure(tr); clf
    figure(1); clf
        set(gcf,'Color','w')
        set(gcf,'Position',[675 135 1085 835])
        for i=1:2
            % Analog Channel MEG data
            t1 = data.time{tr};
            z1 = data.trial{tr}(indch(i),:);
            % Eye Tracking data
            t2 = eyedata(tr).samples(:,1)/1000;
            z2 = eyedata(tr).samples(:,i+1);
            % Interpolate the Eye Tracking data into the MEG timepoints
            z2i = interp1(t2,z2,t1);
            
            fixs_bgn = (eyedata(tr).fixs(:,1) - eyedata(tr).samples0)/1000;
            fixs_end = (eyedata(tr).fixs(:,2) - eyedata(tr).samples0)/1000;
            fixs_bgn_c = (eyedata(tr).fixs(:,1) - eyedata(tr).samples0 - mean(xc_xy(tr,:)))/1000;
            fixs_end_c = (eyedata(tr).fixs(:,2) - eyedata(tr).samples0 - mean(xc_xy(tr,:)))/1000;
            

            subplot(2,2,2*(i-1)+1)
                hold on
                    plot(t2,z2,'k-')
%                     plot(t1,z2i,'r-')
                    
                    YLIMI = ylim;
                    for j=1:size(eyedata(tr).fixs,1)
                        plot([1 1]*fixs_bgn(j), ...
                                YLIMI,'r-')
                        plot([1 1]*fixs_end(j), ...
                                YLIMI,'b-')
                    end
                hold off
                set(gca,'XLim',[min([t1(1) t2(1)]) max([t1(end) t2(end)])])
                ax1 = gca; % current axes
                ylabel(ax1,label_xy{i})
                xlabel(ax1,'ET Time (s)')

            subplot(2,2,2*(i-1)+2)
                hold on
                    plot(t1,z1,'k-')
                    
                    YLIMI = ylim;
                    for j=1:size(eyedata(tr).fixs,1)
                        plot([1 1]*fixs_bgn(j), ...
                                YLIMI,'r-')
                        plot([1 1]*fixs_end(j), ...
                                YLIMI,'b-')
                        plot([1 1]*fixs_bgn_c(j), ...
                                YLIMI,'r--')
                        plot([1 1]*fixs_end_c(j), ...
                                YLIMI,'b--')
                    end
                hold off
                set(gca,'XLim',[min([t1(1) t2(1)]) max([t1(end) t2(end)])])
                ax2 = gca;
                ylabel(ax2,'Amplitude')
                xlabel(ax2,'MEG Time (s)')

        end

    eyedata(tr).megevts.fixOnset     = eyedata(tr).fixs(:,1) - eyedata(tr).samples0 - mean(xc_xy(tr,:));
    eyedata(tr).megevts.fixDuration  = eyedata(tr).fixs(:,3);
    eyedata(tr).megevts.fixPosition  = eyedata(tr).fixs(:,4:5);
    eyedata(tr).megevts.fixRank      = 1:length(eyedata(tr).megevts.fixOnset);

    eyedata(tr).megevts.prevfixDuration   = nan(eyedata(tr).Nfix,1);
    eyedata(tr).megevts.prevSaccDuration  = nan(eyedata(tr).Nfix,1);
    eyedata(tr).megevts.prevSaccAmplitude = nan(eyedata(tr).Nfix,1);
    eyedata(tr).megevts.prevSaccDirection = nan(eyedata(tr).Nfix,1);
    for j=2:eyedata(tr).Nfix
        eyedata(tr).megevts.prevfixDuration(j)   = eyedata(tr).fixs(j-1,3);

        indsac = find(eyedata(tr).sacs(:,2)<eyedata(tr).fixs(j,1),1,'last');
        eyedata(tr).megevts.prevSaccDuration(j)  = eyedata(tr).sacs(indsac,3);

        xini = eyedata(tr).sacs(indsac,4);
        yini = eyedata(tr).sacs(indsac,5);
        xend = eyedata(tr).sacs(indsac,6);
        yend = eyedata(tr).sacs(indsac,7);
        eyedata(tr).megevts.prevSaccAmplitude(j) = sqrt( (xini - xend).^2 + (yini - yend).^2 );
        eyedata(tr).megevts.prevSaccDirection(j) = fun_dirsacc_theta(xini,yini,xend,yend);
    end

    istarget = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,eyedata(tr).info.filenames_target));    
    if strfind(eyedata(tr).info.filenames_target,'object')
        isrelev = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'object'));    
        isirrel = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'face'));        
    elseif strfind(eyedata(tr).info.filenames_target,'face')
        isrelev = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'face'));    
        isirrel = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'object'));        
    end
    isobject    = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'object'));    
    isface      = cellfun(@(x) ~isempty(x), strfind(eyedata(tr).info.filenames_grid,'face'));        

    eyedata(tr).megevts.fixTrialType  = repmat(find(strcmp(eyedata(tr).info.trialtype,{'VS','EX','ME'})),eyedata(tr).Nfix,1);
    eyedata(tr).megevts.fixStimType   = nan(eyedata(tr).Nfix,1);
    eyedata(tr).megevts.fixType   = nan(eyedata(tr).Nfix,1);
    for j=1:eyedata(tr).Nfix
        if ( eyedata(tr).stposrel(j)==0 )
            eyedata(tr).megevts.fixType(j) = 0;
        else
            if isirrel(eyedata(tr).stposrel(j))
                eyedata(tr).megevts.fixType(j) = 1;
            elseif isrelev(eyedata(tr).stposrel(j))
                eyedata(tr).megevts.fixType(j) = 2;
            end
            if istarget(eyedata(tr).stposrel(j))
                eyedata(tr).megevts.fixType(j) = 3;
            end

            if isobject(eyedata(tr).stposrel(j))
                eyedata(tr).megevts.fixStimType(j) = 2;
            elseif isface(eyedata(tr).stposrel(j))
                eyedata(tr).megevts.fixStimType(j) = 1;
            end
        end
    end
    pause(0.1)
end

%% Re-epoch data
% data = 
%          label: {374x1 cell}
%      trialinfo: [160x1 double]
%     sampleinfo: [160x2 double]
%          trial: {1x160 cell}
%           time: {1x160 cell}
%            cfg: [1x1 struct]
%        fsample: 600
%
% In data doesn't seem to be an event variable to chance and re-epoch the
% data. But, on the other hand, I can do it without changing much.
epochLimits_ms  = [-200 250]; % In ms,
epochLimits_smp = round(data.fsample*epochLimits_ms/1000); % In ms,

thrDuration = [250 1000]; % In ms, [200 NaN] means only larger than 200ms
k = 0;
sampleinfo = nan(1000,2);
fixcount = nan(length(data.time),2);
clear trial time megevts
for tr = 1:length(data.time);
    fixcount(tr,1) = eyedata(tr).Nfix;
    if ( eyedata(tr).Nfix > 0 )
        if ~any(isnan(thrDuration))
            ind = find(eyedata(tr).megevts.fixType~=0 & eyedata(tr).megevts.fixDuration >= thrDuration(1) & eyedata(tr).megevts.fixDuration <= thrDuration(2));
        elseif isnan(thrDuration(1))
            ind = find(eyedata(tr).megevts.fixType~=0 & eyedata(tr).megevts.fixDuration <= thrDuration(2));
        elseif isnan(thrDuration(2))
            ind = find(eyedata(tr).megevts.fixType~=0 & eyedata(tr).megevts.fixDuration >= thrDuration(1));
        else
            ind = find(eyedata(tr).megevts.fixType~=0);
        end

        fixcount(tr,2) = sum(~isnan(eyedata(tr).megevts.fixOnset));
        for i = 1:length(ind)
            if ~isnan(eyedata(tr).megevts.fixOnset(ind(i)));
                indt0 = find( data.time{tr} >= eyedata(tr).megevts.fixOnset(ind(i))/1000,1,'first');
                indt = (indt0+epochLimits_smp(1)):(indt0+epochLimits_smp(2));
                
                if (max(indt) <= length(data.time{tr}))
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
k
sampleinfo = sampleinfo(sum(isnan(sampleinfo)')'==0,:);

%%
dataMEG = data;
dataMEG.time       = time;
dataMEG.trial      = trial;
dataMEG.sampleinfo = sampleinfo;
dataMEG.trialinfo  = ones(size(time))'; % My variable of interest
dataMEG.cfg.trials = ones(size(time))'; % My variable of interest

load('labels_CTF269Nott','labels');
scfg.channel = labels;
dataMEG = ft_selectdata(scfg,dataMEG);

bscfg = [];
bscfg.demean            = 'yes';        % 'no' or 'yes', whether to apply baseline correction (default = 'no')
bscfg.baselinewindow    = [-0.2 -0.1];  % [begin end] in seconds, the default is the complete trial (default = 'all')
cfg.reref               = 'yes';        % 'no' or 'yes' (default = 'no')
cfg.refchannel          = 'all';        % cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
cfg.refmethod           = 'avg';        % 'avg' or 'median' (default = 'avg')

[dataMEG] = ft_preprocessing(bscfg, dataMEG);

%%
clear timelock;
for i=1:2
    tmlkcfg = [];
    tmlkcfg.trials = ([megevts.fixStimType] == i & [megevts.fixType]~=3);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'MEG';
    tmlkcfg.lpfilter    = 'yes';
    tmlkcfg.lpfreq      = 45;
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'yes';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataMEG); 
end    

%%
figure(10); clf
    set(gcf,'Position',[80 1 800 970],'Color','w')
    for i=1:2
        subplot(3,4,4*(i-1)+[1:2]);
            hold on
                imagesc(timelock(i).time,1:length(timelock(i).label),timelock(i).avg,[-1.5 1.5]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            set(gca,'XTickLabel',[])
            ylabel('Sensors')
            
        subplot(3,4,4*(i-1)+3);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.5 1.5]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            ft_topoplotER(tpcfg,timelock(i)); %colorbar
            
        subplot(3,4,4*(i-1)+4);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.5 1.5]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            ft_topoplotER(tpcfg,timelock(i)); %colorbar
    end
        subplot(3,4,9:10);
            hold on
                imagesc(timelock(1).time,1:length(timelock(1).label),timelock(1).avg-timelock(2).avg,[-0.5 0.5]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            ylabel('Sensors')
            xlabel('Time (s)')

            
%%
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to 
                                 % evaluate the effect at the sample level
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that 
                                 % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the 
                                 % permutation distribution. 
cfg.minnbchan = 2;               % minimum number of neighborhood channels that is 
                                 % required for a selected sample to be included 
                                 % in the clustering algorithm (default=0).
% cfg.neighbours = neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 100;      % number of draws from the permutation distribution

design = zeros(1,size(timelock(1).trial,1) + size(timelock(2).trial,1));
design(1,1:size(timelock(1).trial,1)) = 1;
design(1,(size(timelock(1).trial,1)+1):(size(timelock(1).trial,1) + size(timelock(2).trial,1)))= 2;

cfg.design = design;             % design matrix
cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)

%%
cfg_neighb        = [];
cfg_neighb.method = 'distance';         
cfg.layout        = 'CTF269Nott.lay';
cfg.channel       = 'all';
neighbours        = ft_prepare_neighbours(cfg_neighb, dataMEG);

cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with 
                                 % which other sensors it can form clusters
ft_neighbourplot(cfg, data);

cfg.channel       = {'MEG'};     % cell-array with selected channel labels
cfg.latency       = [-0.2 0.25]; % time interval over which the experimental 
                                 % conditions must be compared (in seconds)

[stat] = ft_timelockstatistics(cfg, timelock(1), timelock(2));         
                                 