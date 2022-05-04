% Saccade detection parameters
% 
% Please read Engbert & Mergenthaler (2006) for details about the following settings.
% 
% Velocity threshold multiplier
% Threshold multiplier (λ) for saccade detection. The velocity threshold for each component (x,y) of the eye track is set at λ multiples of a median-based estimator for the standard deviation of all samples in the data epoch. The optimal value for λ depends on your measurement setup and eye tracking hardware. For microsaccade detection, values between λ=4 and λ=6 have been used in the literature.
% 
% Minimum saccade duration
% Minimum duration of saccades in sampling points. Movements lasting less than MINDUR samples will not be detected as saccades.
% 
% Smooth raw data to suppress noise?
% If the checkbox is set, eye position data is smoothed before saccade detection (recommended for ET with high sampling rates).
% 
% Compute velocity thresholds globally (across epochs)?
% If the checkbox is set, the same detection thresholds are applied to all data epochs. Velocity thresholds will not be computed individually for each epoch, but globally across all data epochs. If the dataset is still continuous, this setting has no effect (since there is only one "epoch" anyway). Note: This option has been added for the current toolbox. It is still experimental.
% 
% For clusters separated by less than... do the following
% Without additional processing, parts of the same eye movement trajectory may be detected as more than one saccade. For example, there is often an overshoot component seen at the end of saccades ("glissades", cf. Nyström & Holmqvist, BRM, 2010). If eye velocity drops temporarily below the threshold, the glissade would be detected as a new saccade event, even though it is probably better classified as part of the same movement.
% 
% To cluster saccades, enter a value (in milliseconds) in the input field. This value defines the minimum plausible interval (fixation duration) between successive saccades (e.g., 50 ms). The dropdown menu then offers four options on how to treat clusters of saccades that are closer together than this value:
% 
% 1. keep all saccades: Do nothing. Keep all saccades (even if they occur in temporal proximity)
% 2. keep first saccade: Of each temporal cluster of saccades, keep only the first saccade
% 3. keep largest saccade: Of each temporal cluster of saccades, keep only the saccade with the largest distance
% 4. combine into one saccade: Combine all movements of a cluster into one saccade
% Note: Saccade clustering has been added for the current toolbox. Important: Option 4 (combine into one saccade) is still experimental. While be believe it is useful, it has not been thoroughly evaluated, especially for the special case of microsaccade detection! Use it carefully.

for tr=1:length(tr); eyedata(tr).srate   = 500; end
vfac            = 4;
mindur          = round(0.05*srate);
% MEG screen: 385 x 230
% MEG screen: 1920x 1080
degperpixel     = mean([1920/385 1080/230]); 
smooth          = 1;
globalthresh    = 0;
clusterdist     = 50;
clustermode     = 3;
plotfig         = 1;
writesac        = 1;
writefix        = 1;

% [saccade,fixation] = fun_detecteyemovements(eyedata(1).samples,eyedata(1).srate,vfac,mindur,degperpixel,smooth,globalthresh,clusterdist,clustermode,plotfig,writesac,writefix);
