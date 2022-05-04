clear all
close all
clc

cd('~/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_analysis/')
% addpath('/home/juank/toolbox/fieldtrip-20170827')
% ft_defaults;
addpath('my_functions/')
addpath('markus/')
experimentdir       = '/media/backup/Experimentos/Project_EyeMov_EEG_Completely_FreeV/2018_Notts/Data/MEG_recordings';
% experimentdir       = '/home/juank/Experimentos/2018_Notts/MEG_recordings';
directories         = {'/13229-001'};
sessionfilenames    = {'13229-001_MarkusBauer_20180319'};

isubject = 1;
cfg.directory = [experimentdir,directories{isubject},filesep];

cfg.myfname = ['dataEpoched300_',sessionfilenames{isubject},'_merged.mat'];
load([cfg.directory,'tmp/',cfg.myfname],'dataEpoched','recfg')

eyepath = '/home/juank/Dropbox/fERPs_Wheres_Wally/2018_ExpeNotts/MEG_recordings/data/';
subjid = '13229_001';
load([eyepath subjid '/' subjid '_eyedata.mat']);

bscfg = [];
bscfg.demean            = 'yes';        % 'no' or 'yes', whether to apply baseline correction (default = 'no')
bscfg.baselinewindow    = [-0.2 -0.1];  % [begin end] in seconds, the default is the complete trial (default = 'all')
[dataEpoched] = ft_preprocessing(bscfg, dataEpoched);

%% Faces vs Objects analysis: Averaging
clear timelock;
for i=1:2
    tmlkcfg = [];
    tmlkcfg.trials      = ([recfg.megevts.fixStimType] == i & [recfg.megevts.fixType]~=3);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'all';
%     tmlkcfg.lpfilter    = 'no';
    tmlkcfg.lpfilter    = 'yes';
    tmlkcfg.lpfreq      = 45;
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'no';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataEpoched); 
end  

%% Faces vs Objects analysis: Plot fERFs imagesc and topographies
% lycfg = [];
% lycfg.layout = 'CTF269Nott.lay';   
% ft_layoutplot(lycfg, data)
label_cond = {'Faces','Objects'};
figure(10); clf
    set(gcf,'Position',[80 1 800 970],'Color','w')
    for i=1:2
        subplot(3,4,4*(i-1)+[1:2]);
            hold on
                title(label_cond{i})
                imagesc(timelock(i).time,1:length(timelock(i).label),timelock(i).avg,[-1.0 1.0]*10^-13)
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
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
        subplot(3,4,4*(i-1)+4);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
    end
        dtimelock = timelock(1);
        dtimelock.avg = timelock(1).avg-timelock(2).avg;
        subplot(3,4,9:10);
            hold on
                title('Faces-Objects')
                imagesc(timelock(1).time,1:length(timelock(1).label),dtimelock.avg,[-1.0 1.0]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            ylabel('Sensors')
            xlabel('Time (s)')

        subplot(3,4,11);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
        subplot(3,4,12);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
            
%% Faces vs Objects analysis: Plot fERFs waveforms          
label_cond = {'Faces','Objects'};
col_cond    = {[1.0 0.0 0.0],[0.0 0.0 1.0]};
shadow_cond = {[1.0 0.5 0.5],[0.5 0.5 1.0]};

indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLO')),dataEpoched.label)), ...
            find(cellfun(@(x) ~isempty(strfind(x,'MRO')),dataEpoched.label))}; 
label_chan = {'Left Occ','Right Occ'};
% indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLF1')),data.label)), ...
%             find(cellfun(@(x) ~isempty(strfind(x,'MRF1')),data.label))}; 
% label_chan = {'Left Frontal','Right Frontal'};
figure(11); clf
    set(gcf,'Position',[80 1 800 400],'Color','w')
    
    YLIMI = [-10^-13 10^-13];
    for j=1:2
        subplot(1,2,j);
            hold on
                title(label_chan{j});
%                 plot(timelock(i).time,mean(timelock(1).avg(indch{j},:)) - mean(timelock(2).avg(indch{j},:)),...
%                     'Color','k','LineWidth',2)
                hp = nan(2,1);
                for i=1:2
                    [m,e]=fun_meancurve(timelock(i),indch{j});
                    errorPatch = niceBars2(timelock(i).time,m,e,shadow_cond{i},0.5);hold on;
                    hp(i)=plot(timelock(i).time,m,'Color',col_cond{i},'LineWidth',3);
                    set(errorPatch,'EdgeAlpha',0)
                    
                end
                plot([0 0],YLIMI,'k--')
                plot([min(timelock(i).time) max(timelock(i).time)],[0 0],'k--')
            hold off
            set(gca,'XLim',[min(timelock(i).time) max(timelock(i).time)])
            set(gca,'YLim',[YLIMI])
            xlabel('Time (s)')
            ylabel('Amplitude (fT)')
            if (j==2); legend(hp,label_cond,'FontSize',10); legend boxoff; end
    end
   
%% Irrelevant vs Relevant vs Target (VS, {'VS','EX','ME'}) analysis: Averaging
clear timelock;
for i=1:3
    tmlkcfg = [];
    tmlkcfg.trials      = ([recfg.megevts.fixType] == i & [recfg.megevts.fixTrialType]==1);
    tmlkcfg.vartrllength= 0;
    tmlkcfg.channel     = 'all';
%     tmlkcfg.lpfilter    = 'no';
    tmlkcfg.lpfilter    = 'yes';
    tmlkcfg.lpfreq      = 45;
    tmlkcfg.keeptrials  = 'yes';
    tmlkcfg.covariance  = 'no';
    timelock(i) = ft_timelockanalysis(tmlkcfg,dataEpoched); 
end  

%% Irrelevant vs Relevant vs Target (VS, {'VS','EX','ME'}) analysis: Plot fERFs imagesc and topographies
% lycfg = [];
% lycfg.layout = 'CTF269Nott.lay';   
% ft_layoutplot(lycfg, data)
label_cond = {'Irrelev','Relev','Targets'};
indsel = [2 3];
figure(30); clf
    set(gcf,'Position',[80 1 800 970],'Color','w')
    for ii=1:2
        i = indsel(ii);
        subplot(3,4,4*(ii-1)+[1:2]);
            hold on
                title(label_cond{i})
                imagesc(timelock(i).time,1:length(timelock(i).label),timelock(i).avg,[-1.0 1.0]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            set(gca,'XTickLabel',[])
            ylabel('Sensors')
            
        subplot(3,4,4*(ii-1)+3);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
        subplot(3,4,4*(ii-1)+4);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,timelock(i)); %colorbar
            hold off
    end
        dtimelock = timelock(indsel(1));
        dtimelock.avg = timelock(indsel(1)).avg-timelock(indsel(2)).avg;
        subplot(3,4,9:10);
            hold on
                title([label_cond{indsel(1)} '-' label_cond{indsel(2)}])
                imagesc(timelock(1).time,1:length(timelock(1).label),dtimelock.avg,[-1.0 1.0]*10^-13)
                plot([0 0],[1 length(timelock(i).label)],'k--')
            hold off
            hc = colorbar;
            ylabel(hc,'Amplitude (fT)')
            set(gca,'XLim',[timelock(i).time(1) timelock(i).time(end)])
            set(gca,'YLim',[0.5 length(timelock(i).label)+0.5])
            ylabel('Sensors')
            xlabel('Time (s)')

        subplot(3,4,11);
            tpcfg = [];                            
            tpcfg.xlim = [0.080 0.120];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'avg'; % the default 'avg' is not present in the data
            hold on
                title('80-120ms')
                
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
        subplot(3,4,12);
            tpcfg = [];                            
            tpcfg.xlim = [0.15 0.20];                
            tpcfg.zlim = [-1.0 1.0]*10^-13;                
            tpcfg.marker = 'off';
            tpcfg.comment = 'no'; 
            tpcfg.layout = 'CTF269Nott.lay';            
%             tpcfg.parameter = 'individual'; % the default 'avg' is not present in the data
            hold on
                title('150-200ms')
                ft_topoplotER(tpcfg,dtimelock); %colorbar
            hold off
            
%% Irrelevant vs Relevant vs Target (VS, {'VS','EX','ME'}) analysis: Plot fERFs waveforms          
label_cond = {'Irrelev','Relev','Targets'};
col_cond    = {[1.0 0.0 0.0],[0.0 1.0 0.0],[0.0 0.0 1.0]};
shadow_cond = {[1.0 0.5 0.5],[0.5 1.0 0.5],[0.5 0.5 1.0]};

indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLO')),dataEpoched.label)), ...
            find(cellfun(@(x) ~isempty(strfind(x,'MRO')),dataEpoched.label))}; 
label_chan = {'Left Occ','Right Occ'};
% indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLF')),dataEpoched.label)), ...
%             find(cellfun(@(x) ~isempty(strfind(x,'MRF')),dataEpoched.label))}; 
% label_chan = {'Left Frontal','Right Frontal'};
figure(11); clf
    set(gcf,'Position',[80 1 800 400],'Color','w')
    
    YLIMI = [-2*10^-13 2*10^-13];
    for j=1:2
        subplot(1,2,j);
            hold on
                title(label_chan{j});
%                 plot(timelock(i).time,mean(timelock(1).avg(indch{j},:)) - mean(timelock(2).avg(indch{j},:)),...
%                     'Color','k','LineWidth',2)
                hp = nan(3,1);
                for i=1:3
                    [m,e]=fun_meancurve(timelock(i),indch{j});
                    errorPatch = niceBars2(timelock(i).time,m,e,shadow_cond{i},0.5);hold on;
                    hp(i)=plot(timelock(i).time,m,'Color',col_cond{i},'LineWidth',3);
                    set(errorPatch,'EdgeAlpha',0)
                    
                end
                plot([0 0],YLIMI,'k--')
                plot([min(timelock(i).time) max(timelock(i).time)],[0 0],'k--')
            hold off
            set(gca,'XLim',[min(timelock(i).time) max(timelock(i).time)])
            set(gca,'YLim',[YLIMI])
            xlabel('Time (s)')
            ylabel('Amplitude (fT)')
            if (j==2); legend(hp,label_cond,'FontSize',10); legend boxoff; end
    end
   
%% Irrelevant vs Relevant vs Target (VS, {'VS','EX','ME'}) analysis: Plot fERFs waveforms          
label_cond = {'Irrelev','Relev','Targets'};
col_cond    = {[1.0 0.0 0.0],[0.0 1.0 0.0],[0.0 0.0 1.0]};
shadow_cond = {[1.0 0.5 0.5],[0.5 1.0 0.5],[0.5 0.5 1.0]};

indch = {   find(cellfun(@(x) ~isempty(strfind(x,'MLF')),dataEpoched.label)), find(cellfun(@(x) ~isempty(strfind(x,'MRF')),dataEpoched.label)); ... 
            find(cellfun(@(x) ~isempty(strfind(x,'MLC')),dataEpoched.label)), find(cellfun(@(x) ~isempty(strfind(x,'MRC')),dataEpoched.label)); ... 
            find(cellfun(@(x) ~isempty(strfind(x,'MLP')),dataEpoched.label)), find(cellfun(@(x) ~isempty(strfind(x,'MRP')),dataEpoched.label)); ... 
            find(cellfun(@(x) ~isempty(strfind(x,'MLO')),dataEpoched.label)), find(cellfun(@(x) ~isempty(strfind(x,'MRO')),dataEpoched.label))}; 
label_chan = {'FL','FR'; 'CL', 'CR'; 'PL', 'PR'; 'OL', 'OR'};
YLIMI = {   [-10^-13 10^-13],[-10^-13 10^-13];[-10^-13 10^-13],[-10^-13 10^-13];...
            [-10^-13 10^-13],[-10^-13 10^-13];2*[-10^-13 10^-13],2*[-10^-13 10^-13]};
figure(11); clf
    set(gcf,'Position',[80 1 800 900],'Color','w')
    
    for j=1:4
        for k=1:2
            subplot(4,2,2*(j-1)+k);
                hold on
                    title(label_chan{j,k});
                    hp = nan(3,1);
                    for i=3:-1:1
                        [m,e]=fun_meancurve(timelock(i),indch{j,k});
                        errorPatch = niceBars2(timelock(i).time,m,e,shadow_cond{i},0.5);hold on;
                        hp(i)=plot(timelock(i).time,m,'Color',col_cond{i},'LineWidth',3);
                        set(errorPatch,'EdgeAlpha',0)

                    end
                    plot([0 0],YLIMI{j,k},'k--')
                    plot([min(timelock(i).time) max(timelock(i).time)],[0 0],'k--')
                hold off
                set(gca,'XLim',[min(timelock(i).time) max(timelock(i).time)])
                set(gca,'YLim',[YLIMI{j,k}])
                xlabel('Time (s)')
                ylabel('Amplitude (fT)')
                if (2*(j-1)+k==8); legend(hp,label_cond,'FontSize',10); legend boxoff; end
        end
    end
   