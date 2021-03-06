% (2018-03-26, JK) Adapted from Dimigen's EYE-EEG

% detecteyemovements() - detect saccades & fixations in eye tracking data. 
%                Saccade detection is based on the algorithm by 
%                Engbert & Mergenthaler (2006). Saccades are defined as 
%                (monocular or binocular) outliers in 2D velocity space.
%                Velocity thresholds for saccade detection are determined 
%                adaptively as a multiple of the (median-based) SD of all 
%                data samples in the epoch. Fixations are defined as 
%                intervals in-between saccades. Eye movements can be added 
%                as new events to EEGLAB's event structure. For various 
%                other options, see below.
%
% Usage:
%   >> EEG = detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,...
%            degperpixel,smooth,globalthresh,clusterdist,clustermode,
%            plotfig,writesac,writefix)
%
% Required inputs:
%   EEG          - [string] EEG struct, also containing synchronized eye 
%                  tracking data (see pop_importeyetracker)
%   left_eye_xy  - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of left eye. Leave empty []
%                  if the left eye was not recorded.
%   right_eye_xy - [vector of two channel indices], 
%                  specifying channel indices of X- (first value) and 
%                  Y-component (second value) of right eye. Leave empty []
%                  if the right eye was not recorded.
%   vfac         - [double] velocity factor ("lambda") to determine 
%                  the velocity threshold for saccade detection 
%                  (cf. Engbert & Mergenthaler, 2006)
%   mindur       - [integer] minimum saccade duration (in samples)
%                  (cf. Engbert & Mergenthaler, 2006)
%   degperpixel  - [double] visual angle of one screen pixel
%                  if this value is left empty [], saccade characteristics 
%                  are reported in the original data metric (pixel?) 
%                  instead of in degrees of visual angle
%   smooth       - [0/1] if set to 1, the raw data is smoothed over a 
%                  5-sample window to suppress noise
%                  noise. Recommended for high native ET sampling rates.
%   globalthresh - [0/1]. Use the same thresholds for all epochs? 
%                  0: Adaptive velocity thresholds are computed 
%                  individually for each data epoch. 
%                  1: Adaptive velocity thresholds are first computed for 
%                  each epoch, but then the mean thresholds are applied to 
%                  each epochs (i.e. same detection parameters are used for
%                  all epochs). Setting is irrelevant if the input data is 
%                  still continuous (= only one data epoch).
%   clusterdist  - [integer] value in sampling points that defines the
%                  minimum allowed fixation duration between two saccades. 
%                  If the off- and onsets of two temp. adjacent sacc. are 
%                  closer together than 'clusterdist' samples, these 
%                  saccades are regarded as a "cluster" and treated 
%                  according to the 'clustermode' setting (see below).
%                  clusterdist is irrelevant if clustermode == 1.
%   clustermode  - [1,2,3,4]. Integer between 1 and 4.
%                  1: keep all saccades, do nothing
%                  2: keep only first saccade of each cluster
%                  3: keep only largest sacc. of each cluster
%                  4: combine all movements into one (longer) saccade
%                     this new saccade is defined as the movement that 
%                     occurs between the onset of the 1st saccade in the
%                     cluster and the offset of the last sacc. in cluster
%                     WARNING: CLUSTERMODE 4 is experimental and untested!  
%   plotfig      - [0/1] Show a figure with eye movement properties?
%                  0: do not plot a figure. 
%                  1: plot a figure displaying properties of detected 
%                  saccades & fixations
%   writesac     - [0/1]: Add saccades to EEG.event?
%                  0: detect saccades, but do not store them in EEG.event. 
%                  1: add detected saccades as new events to EEG.event.
%   writefix     - [0/1]: Add fixations to EEG.event?
%                  0: detect fixations, but do not add them to EEG.event.
%                  1: add detected fixations as new events to EEG.event.
%                  Note: It is recommended to first test the parameters of
%                  saccade/fixation detection without adding events.
%                  For this, set writesac and writefix to 0.
%
% Outputs:
%   EEG         - EEG structure. If writesac or writefix were set to 1, 
%                 the EEG structure (EEG.event/EEG.urevent/EEG.epoch)
%                 will contain additional "saccade" and "fixation" events
%                 with their respective properties
%
% See also: vecvel, velthresh, microsacc_plugin, binsacc, saccpar, 
%           mergesacc, addevents
%
%
% An example call of the function might look like this: 
% >> EEG = detecteyemovements(EEG,[],[33 34],6,4,0.037,1,0,25,4,1,1,0)
%
% In this example, the eye position data for the right eye is stored in 
% channels 33 (horiz.) and 34 (vertical). The left eye was not recorded. 
% The velocity threshold is set to 6 times the (median-based) 
% SD of all velocity samples in the epoch. The minimum duration of
% saccades to be detected is 4 samples. In the experiment, one screen 
% pixel corresponded to 0.037 degrees of visual angle. 
% The raw data is smoothed prior to saccade detection (smooth: 1). 
% Adaptive velocity thresholds (X and Y-threshold for each eye) are 
% determined individually for each data epoch (globalthresh: 0). For saccades 
% separated by fixations of less than 25 samples, only the first saccade 
% is kept (clusterdist: 25, clustermode: 2). A figure with the 
% results is plotted. Detected saccades are stored as new events in 
% EEG.event, but fixations are not stored.
% 
% The eye movement detection is based on:
%
% Engbert, R., & Kliegl, R. (2003). Microsaccades uncover the orientation
% of covert attention. Vision Research, Vol. 43, 1035-1045
%
% ...as well as...
%
% Engbert, R., & Mergenthaler, K. (2006). Microsaccades are triggered by 
% low retinal image slip, PNAS, Vol. 103 (18), 7192-7197
%
% Author: od
% Copyright (C) 2009-2017 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de 

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

% function EEG = fun_detecteyemovements(EEG,left_eye_xy,right_eye_xy,vfac,mindur,degperpixel,smooth,globalthresh,clusterdist,clustermode,plotfig,writesac,writefix)
function [saccade,fixation] = fun_detecteyemovements(samples,srate,vfac,mindur,degperpixel,smooth,globalthresh,clusterdist,clustermode,plotfig,writesac,writefix)
allsac = [];
allfix = [];

nsample     = size(samples,1);
nbadsmp     = 0;

% warn message if back-to-back saccades are detected
clusterwarning = false;

% critical bugfix 2013-10-01, by OD
% due to bug in third-party function "smoothdata":
% function smoothdata() was removed from the toolbox
% function vecvel() was updated to incorporate different levels of smoothing:
% options:
% - smoothlevel 0: no smoothing, simple diff()
% - smoothlevel 1: 3-point window
% - smoothlevel 2: 5-point window
if smooth
    smoothlevel = 2; % 5-point smoothing
else
    smoothlevel = 0; % no smoothing
end

%% screen feedback
% summary of detection parameters
fprintf('\n--------------------------------------------------------------------')
fprintf('\nDetecting saccades after Engbert & Mergenthaler (2006)\n')
fprintf('\nVelocity threshold factor (vfac):  %.2f SD',vfac);
fprintf('\nMinimum saccade duration (mindur): %.2f samples (%.2f ms)',mindur,mindur*1000/srate);
if ~isempty(degperpixel) | isnan(degperpixel) % bugfix 2016-11-12 by OD: added case if degperpixel = NaN (from GUI input) 
    fprintf('\nVisual angle per screen pixel:     %f???',degperpixel);
    metric = 'deg';
else
    fprintf('\nWARNING: No input provided for degperpixel!\nSpatial saccade properties are given in original metric (pixel?)');
    degperpixel = 1;
    metric = 'pix';
end
if smooth, fprintf('\n-- Raw data is smoothed in 5-sample window.'); else fprintf('\n-- Raw data is not smoothed.'); end
fprintf('\n-- Treatment of saccade clusters:');
switch clustermode    
    case 1
        fprintf('\n\tAll saccades are kept');
    case 2
        fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/srate));
        fprintf('\n\tFirst sacc. of each cluster is kept.');
    case 3
        fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/srate));
        fprintf('\n\tLargest sacc. of each cluster is kept.');
    case 4
        fprintf('\n\tSaccades separated by fixations < %i ms are clustered.',clusterdist*(1000/srate));
        fprintf('\n\tClusters are combined into one saccade.');
    otherwise
        error('%s(): Unknown input for clustermode, should be: 1,2,3,4.',mfilename)
end
if plotfig, fprintf('\n-- A figure with eye movement properties is plotted.'); end        
fprintf('\n')

%% pre-compute saccade velocity thresholds for all epochs
z = samples(:,[2 3]);
vz = vecvel(z,srate,smoothlevel);
[z_msdx, z_msdy] = velthresh(vz);

%% detect saccades & fixations
sac = [];
z = samples(:,[2 3]);
% bad/missing samples in eye track?
badsmp = sum(sum(z<=0)); if badsmp > 0, badepochs = 1; nbadsmp = nbadsmp + badsmp; end
vz = vecvel(z,srate,smoothlevel); % get eye velocities
% detect monocular saccades
if globalthresh % use precomputed velocity thresholds (mean of all epochs)
    sacZ = microsacc_plugin(z,vz,vfac,mindur,mean(z_msdx),mean(z_msdy));
else % compute velocity thresholds from this epoch only
    sacZ = microsacc_plugin(z,vz,vfac,mindur,z_msdx,z_msdy);        
end

sac = sacZ; clear sacL;
sac = saccpar([sac sac]);
sac = mergesacc(sac,z,clusterdist,clustermode);

%% update various saccade metrics


if ~isempty(sac)
    % define saccade duration as difference between saccade offset and
    % saccade onset sample. In saccpar(), monocular saccade durations
    % of both eyes are averaged, leading to uneven values (e.g.: 10.5
    % samples) different from the difference between onset and offset
    % values (which are the monocular extremes).
    % Instead: use difference between offset and onset
    sac(:,3) = sac(:,2)-sac(:,1)+1;

    % report saccade velocity/distance/amplitude as visual angles
    sac(:,[5 6 8]) = sac(:,[5 6 8]) .* degperpixel;

    % report saccade angles in degree rather than radians
    sac(:,[7 9]) = sac(:,[7 9]) * 180/pi;

    % add index of corresp. data epoch
%     sac(:,10) = e;

    % store screen location for start and end of saccade
    gazexy = z;

    % get position immediatly before sacc. onset and after sacc. offset
    startsmp = sac(:,1)-1; endsmp = sac(:,2)+1;
    if startsmp(1) < 1, startsmp(1,1) = 1; end
    if endsmp(end) > size(gazexy,2), endsmp(end) = size(gazexy,2); end
    sac(:,11) = gazexy(startsmp,1);
    sac(:,12) = gazexy(startsmp,2);
    sac(:,13) = gazexy(endsmp  ,1);
    sac(:,14) = gazexy(endsmp  ,2);
end

% columns of [sac]:
% 1: saccade onset (sample)
% 2: saccade offset (sample)
% 3: duration (samples)
% 4: delay between eyes (samples)
% 5: vpeak (peak velocity)
% 6: saccade "distance" (eucly. dist. between start and endpoint)
% 7: saccade angle (based on saccade "distance")
% 8: saccade "amplitude" (eucly. dist. of min/max in full saccade trajectory)
% 9: saccade angle (based on saccade "amplitude")
%10: index of corresponding data epoch (1 in case of contin. data)
%11: horizontal (x) gaze position before start of saccade (pixel)
%12: vertial (y) gaze position before start of saccade (pixel)
%13: horizontal (x) gaze position after end of saccade (pixel)
%14: vertial (y) gaze position after end of saccade (pixel)

%% get fixations
nsac = size(sac,1);
fix = [];
if ~isempty(sac)
    if nsac > 1
        for f = 1:nsac-1
            fix(f,1) = sac(f,2)+1;
            fix(f,2) = sac(f+1,1)-1;

            % catch special case: if Engbert algorithms are applied 
            % without any saccade clustering [>> mergesacc()] there can
            % be back-to-back saccades with intervening "fixations" of
            % zero sample duration. Catch this by setting the duration
            % of these fixations to one sample
            if fix(f,1) > fix(f,2)
                fix(f,1) = fix(f,2); % 1-sample fixation
                clusterwarning = true;
            end

        end
    end
    % if epoch does not begin with saccade, add first fixation
    if sac(1,1) > 1
        fix = [[1 sac(1,1)-1]; fix];
    end
    % if epoch does not end with saccade, add last fixation
    if sac(end,2) < nsample
        fix = [fix;[sac(end,2)+1 nsample]];
    end

    for f = 1:size(fix,1)
        % fixation duration (in samples!)
        fix(f,3) = fix(f,2)-fix(f,1)+1;

        % mean fix. position: left eye
        fix(f,4) = mean( z(fix(f,1):fix(f,2),1) );
        fix(f,5) = mean( z(fix(f,1):fix(f,2),2) );
        fix(f,6) = NaN;
        fix(f,7) = NaN;

        % binocular fixation position
        fix(f,8) = nanmean(fix(f,[4 6]));
        fix(f,9) = nanmean(fix(f,[5 7]));
    end
%     fix(:,10) = e; % add index of corresp. data epoch

%     % recompute latencies of eye movement events 
%     % (only necessary for epoched data)
%     offset = (e-1)*nsample;
%     sac(:,[1 2]) = sac(:,[1 2])+offset;
%     % special case: single sacc. lasts entire epoch
%     if ~isempty(fix) 
%         fix(:,[1 2]) = fix(:,[1 2])+offset;
%     end
end

% columns of [fix]:
% 1: fixation onset (sample)
% 2: fixation offset (sample)
% 3: duration (samples)
% 4: mean fix position (L eye X)
% 5: mean fix position (L eye Y)
% 6: mean fix position (R eye X)
% 7: mean fix position (R eye Y)
% 8: mean fix position (L/R average X)
% 9: mean fix position (L/R average Y)
%10: index of corresponding data epoch (1 in case of contin. data)

% slow, but simple:
allsac = [allsac;sac];
allfix = [allfix;fix];

%% remove artificial eye movements caused by boundaries (data breaks)
% (2018-03-26, JK) I've the data already epoched. But I've to think in a
% substitute for this.
%
% % Remove all EMs whose onset is detected in temporal proximity to boundary. 
% % Otherwise, data breaks will likely result in additional fake sacc./fix.
% % Applies only if eye movements are detected in continuous data.
% if nepochs == 1
% 
%     ix_bnd = find(cellfun(@(x) strcmp(x,'boundary'),{EEG.event.type})); % bug fix: now robust against numeric types
%     
%     % any data breaks?
%     % (due to manual rejections or function pop_rej_eyecontin)
%     if ~isempty(ix_bnd)
%         
%         % minimum distance from data break in milliseconds (hard-coded)
%         BOUNDDIST_MS = 50;
%         BOUNDDIST    = round(BOUNDDIST_MS / (1000/EEG.srate)); 
%         % data break latencies
%         bound_lats   = round([EEG.event(ix_bnd).latency]); 
% 
%         % mark all samples close to data break
%         boundvector  = zeros(1,EEG.pnts);
%         for b = 1:length(bound_lats)
%             lowr = bound_lats(b)-BOUNDDIST;
%             uppr = bound_lats(b)+BOUNDDIST;
%             if lowr <= 0, lowr = 1; end
%             if uppr > EEG.pnts, uppr = EEG.pnts; end
%             boundvector(lowr:uppr) = 1;
%         end
%         nearboundsmp = find(boundvector);
%         
%         % option 1: event onset is close to boundary
%         % fakesac = find(ismember(allsac(:,1),nearboundsmp));
%         % fakefix = find(ismember(allfix(:,1),nearboundsmp));
%         
%         % option 2: event on- or offset is close to boundary
%         fakesac = find(ismember(allsac(:,1),nearboundsmp) | ismember(allsac(:,2),nearboundsmp));
%         fakefix = find(ismember(allfix(:,1),nearboundsmp) | ismember(allfix(:,2),nearboundsmp));
% 
%         allsac(fakesac,:) = [];
%         allfix(fakefix,:) = [];
%         fprintf('\n--------------------------------------------------------------------');
%         fprintf('\nFound %i data breaks (boundary events) in the continuous data',length(ix_bnd));
%         fprintf('\nRemoving eye movements that might be artifacts of data breaks:');
%         fprintf('\nRemoved %i saccades  < %i ms away from a boundary',length(fakesac),BOUNDDIST_MS);
%         fprintf('\nRemoved %i fixations < %i ms away from a boundary',length(fakefix),BOUNDDIST_MS);
%         % future versions: ix_bad = find(ismember({EEG.event.type},'badeye'));
%     end
% end
% 
% if clusterwarning && clustermode == 1
%     fprintf('\n--------------------------------------------------------------------');
%     fprintf('\n\n*** WARNING! ***\nDetected pairs of saccades that ended/started on adjacent data samples.');
%     fprintf('\nYou should probably use saccade clustering! (clustermode: 2,3,4)');
% end
% 

%% user feedback: are there "bad" epochs containing values <= 0?
% if sum(badepochs)>0
%     fprintf('\n--------------------------------------------------------------------');
%     fprintf('\n*** WARNING! ***\n%i of %i epochs contained gaze position values <= 0',sum(badepochs),nepochs);
%     fprintf('\n-- Total number of samples <= 0: %i',nbadsmp);
%     fprintf('\n-- Did you reject data with out-of-range values (''Reject data based on eyetrack'')?');
%     fprintf('\n-- Did you subtract a baseline from eye channels (explaining negative values)?');
%     fprintf('\n-- Note: any blinks/bad data will be erroneously detected as saccades/fixations!');
% end

%% user feedback: saccade & fixation detection
fprintf('\n--------------------------------------------------------------------');
fprintf('\nVelocity thresholds used'); %if nepochs > 1, fprintf(' (mean across epochs):'); end;
fprintf('\n\t(Only) eye.  Horiz.: %.2f %s/s. Vert.: %.2f %s/s',mean(z_msdx*vfac*degperpixel),metric,mean(z_msdy*vfac*degperpixel),metric);
fprintf('\n--------------------------------------------------------------------')
if ~isempty(allsac)
    fprintf('\n%i saccades detected:',size(allsac,1));
    fprintf('\n\tMedian amplitude:  %.2f %s',median(allsac(:,6)),metric);
    fprintf('\n\tMedian duration:   %.2f ms',median(allsac(:,3))*1000/srate);
    fprintf('\n\tMedian peak veloc. %.2f %s/s',median(allsac(:,5)),metric);
end
if ~isempty(allfix)
    fprintf('\n%i fixations detected:',size(allfix,1));
    fprintf('\n\tMedian duration:   %.2f ms',median(allfix(:,3))*1000/srate);
    fprintf('\n\tMedian fix. pos. (only)  eye: Horiz.: %.2f px. Vert.: %.2f px',median(allfix(:,4)),median(allfix(:,5)));
end
fprintf('\n')

%% plot figure with eye movements properties
if plotfig
    fprintf('\n--------------------------------------------------------------------')
    fprintf('\nPlotting eye movement properties of left eye...')
    ploteyemovements(allsac(:,6),allsac(:,5),allsac(:,7),allfix(:,3)*1000/srate,allfix(:,4),allfix(:,5),metric);
end

%% write eye movements to EEG.event & EEG.urevent
if writefix || writesac   
    % write saccades
    if writesac
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nAdding %i saccades to EEG.event...\n',size(allsac,1))
        saccade.latency         = allsac(:,1);
        saccade.duration        = allsac(:,3);
        saccade.sac_vmax        = allsac(:,5);
        saccade.sac_amplitude   = allsac(:,6);
        saccade.sac_angle       = allsac(:,7);
        saccade.epoch           = allsac(:,10);
        saccade.sac_startpos_x  = allsac(:,11);
        saccade.sac_startpos_y  = allsac(:,12);
        saccade.sac_endpos_x    = allsac(:,13);
        saccade.sac_endpos_y    = allsac(:,14);
    end
    
    % write fixations
    if writefix
        fprintf('\n--------------------------------------------------------------------')
        fprintf('\nAdding %i fixations to EEG.event...\n',size(allfix,1))
        fixation.latency        = allfix(:,1);
        fixation.duration       = allfix(:,3);
        fixation.fix_avgpos_x   = allfix(:,8);
        fixation.fix_avgpos_y   = allfix(:,9);
%         fixation.epoch          = allfix(:,10);
    end
    fprintf('--------------------------------------------------------------------')
end
fprintf('\nDone.\n\n')
end