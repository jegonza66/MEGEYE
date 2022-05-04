%% Synchronize the MEG time with the ET time
function [eyedata,xc_xy] = fun_estimating_lag(data,eyedata,verb)
if (nargin<2)
    verb=0;
end
indch = find(cellfun(@(x) ~isempty(strfind(x,'UADC')),data.label));
data.label(indch);

%% Estimate max of cross-correlation
xc_xy       = nan(length(data.time),2);
if (verb); label_xy    = {'Horizontal Position (pxls)','Vertical Position (pxls)'}; end
for tr = 1:length(data.time);
    if (verb); 
    %     figure(tr); clf
        figure(1); clf;
            set(gcf,'Color','w','Position',[675 135 1085 835])
    end
    for i=1:2

        if (verb); 
            [xc_xy(tr,i),mxcf,xcf,tlags,t1,z1,t2,z2,z2i]=...
                fun_xcorr_lag(tr,i,data,indch,eyedata);
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
                    plot(xc_xy(tr,i),mxcf,'ro');
                    title(xc_xy(tr,i))
                hold off
                xlabel('lags [ms]'), ylabel('corr')
                xlim([-100 100])
            pause(0.1)
        else
            xc_xy(tr,i) = fun_xcorr_lag(tr,i,data,indch,eyedata);
            
        end
    end
end

%% 
% * Correct samples by a fixed time.
% * Calculate fixation/saccades from samples.
% * Change fixation times to MEG time points.

if (verb); label_xy    = {'Horizontal Position (pxls)','Vertical Position (pxls)'}; end
for tr = 1:length(data.time);
    if (verb);
        %     figure(tr); clf
        figure(1); clf
            set(gcf,'Color','w','Position',[675 135 1085 835])
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
        pause(0.1)
    end

    if ~isempty(eyedata(tr).blinks)
        eyedata(tr).Nblink       = size(eyedata(tr).blinks,1);
        eyedata(tr).megevts.blinkOnset   = eyedata(tr).blinks(:,1) - eyedata(tr).samples0 - mean(xc_xy(tr,:));
        eyedata(tr).megevts.blinkDuration= eyedata(tr).blinks(:,3);
    else
        eyedata(tr).Nblink       = 0;
        eyedata(tr).megevts.blinkOnset   = [];
        eyedata(tr).megevts.blinkDuration= [];
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
end
end
      
                                 