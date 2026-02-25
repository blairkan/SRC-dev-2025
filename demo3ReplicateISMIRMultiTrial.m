% demo3ReplicateISMIRMultiTrial.m
% --------------------------------
% Blair dev version 2/18/2026
% stim2eeg demo
%
% perform a SRC analysis on sample data from a single participant 
% viewing a ~5.5 minute clip from the film "Dog Day Afternoon"
%
% Workflow
% 1. Load input data
% 2. Prepare stim feature and input EEG data
% 3. Create stim feature toeplitz matrix
% 4. Call core function computeCCA.m
% 5. Compute SRC --> moved outside of core function
% 6. Visualize outputs


% History
% 2/18/2026 - adapted from demo2ComputeCCA
% 2/3/2026 - adapted from demo.m

clear all; close all; clc

%% Download ISMIR data

% The following script will download any ISMIR files into the Data/ 
% directory if they are not already there. For any file that is already 
% there, the script will skip downloading a second time. The downloaded
% files are ignored by git. 

downloadISMIRData

%% Load input data

%%% Load stim feature files

% Gang 2017 paper: 'We performed all subsequent analyses using PC1, as 
% well as two individual features. RMS and spectral flux were chosen as
% they reflect amplitude envelope and timbre, respectively, and have been 
% used in previous studies mapping music stimulus features to brain 
% responses.'

% In this file, 'vars' is a struct with named fields for each feature. Each
% feature is a column vector
tempX1 = load('21_Features_80Hz.mat', 'fs', 'vars'); 

% In this file, 'PC' is a column vector feature
tempX2 = load('21_PC1_80Hz.mat', 'fs', 'PC');

%%% Load EEG file

% In this file, 'data21' is a [space x time x trial] matrix
tempY = load('song21_all_Imputed.mat', 'fs', 'data21'); 

%%% Check and extract sampling rates
assert(isequal(tempX1.fs, tempX2.fs), 'Stim feature sampling rates are not equal!')
fsStim = tempX1.fs; % 80
fsEEG = tempY.fs; % 125

%%% Extract stim features of interest
stimRMS = tempX1.vars.RmsEnergy;
stimFlux = tempX1.vars.Flux; 
stimPC1 = tempX2.PC; 

%%% Extract EEG data matrix
eeg = tempY.data21; % [space x time x trial] matrix, no NaNs
[nSpace, nTime, nTrial] = size(eeg); 
disp(['Input data has ' num2str(nSpace) ' channels, ' num2str(nTime) ...
    ' time points, and ' num2str(nTrial) ' trials of data.'])

clear temp*

%% Downsample and DC correct the EEG

%%% Downsample EEG to sampling rate of the stim feature
% General notes:
% - We will generally downample to the lower sampling rate, since with 
%   linear methods there is no value add (i.e., additional information) 
%   from upsampling to the higher sampling rate
% - As is the case here, the stim feature often has the lower sampling rate
% - Therefore, we are downsampling the EEG to the sampling rate of the stim
%   feature

% Iterate through the trials to resample each trial and save the outputs 
% into a new 3D matrix. 
% - Note that EEG input to the resample function should have  channels as 
%   columns, which is why we are transposing the input 
%   [space x time x trial] matrices. 
% - We also transpose the outputs of the resample function call so that the 
%   resulting 3D matrix is still be [space x time x trial].
for i = 1:nTrial
    eegDown(:, :, i) = resample(eeg(:, :, i)', fsStim, fsEEG)';
end

%%% DC Correct the EEG
% Now we will iterate through the trials again and DC correct each trial,
% saving the outputs into a new 3D matrix.
% - the dcCorrect function is one of Blair's functions. It operates on
%   [space x time] input matrices. 

for i = 1:nTrial    
    eegDC(:, :, i) = dcCorrect(eegDown(:, :, i));
end

% Get data size again
[nSpace, nTime, nTrial] = size(eegDC);

%% z-score the stim feature and match its length to the EEG

%%% Normalize the feature because correlation doesn't care about scale
zStimRMS = zscore(stimRMS);
zStimFlux = zscore(stimFlux); 
zStimPC1 = zscore(stimPC1);

%%% Adjust length of stim feature to match length of EEG
% - It's common for the downsampled EEG and the stimulus to have slightly
%   different lengths, so we fix this here
% - Currently the procedure is to trim or zero pad the stim feature input 
%   to match the EEG, and never adjust the EEG. We could adjust this later
%   if needed (to e.g., always trim to the shortest input)

%%% Check each of the features and adjust if needed
% RMS feature
if nTime > length(zStimRMS) % If EEG is longer than the feature
    disp([newline 'RMS: Zero-padding to length of EEG'])
    nZero = nTime - length(zStimRMS);
    zStimRMS = [zStimRMS; zeros(nZero,1)];
elseif length(zStimRMS) > nTime % If feature is longer than EEG
    disp('RMS: Trimming to length of EEG')
    zStimRMS = zStimRMS(1:nTime);
end

% Flux feature
if nTime > length(zStimFlux) % If EEG is longer than the feature
    disp([newline 'Flux: Zero-padding to length of EEG'])
    nZero = nTime - length(zStimFlux);
    zStimFlux = [zStimFlux; zeros(nZero,1)];
elseif length(zStimFlux) > nTime % If feature is longer than EEG
    disp('Flux: Trimming to length of EEG')
    zStimFlux = zStimFlux(1:nTime);
end

% PC1 feature
if nTime > length(zStimPC1) % If EEG is longer than the feature
    disp([newline 'PC1: Zero-padding to length of EEG'])
    nZero = nTime - length(zStimPC1);
    zStimPC1 = [zStimPC1; zeros(nZero,1)];
elseif length(zStimPC1) > nTime % If feature is longer than EEG
    disp('PC1: Trimming to length of EEG')
    zStimPC1 = zStimPC1(1:nTime);
end

%% Create matrix of time-shifted versions of stim feature
% Now we need to create a convolution matrix from the one-dimensional
% feature time series. This allows us to temporally filter the stimulus
% time series using a matrix-vector product.
% - createFeatureToeplitzMatrix() is Blair's function, which converts a 
% vector to a convolution matrix of successive shifts. It functions in the 
% same way as the tplitz() function that was originally called here, but
% avoids out-of-memory errors that can happen with tplitz().

%%% Before calling the function, we need to specify how long we want the
% filter to be.
% - Here we set this to one second, which means we look back up to 1 second 
%   in time. You can think of this as the maximum delay between the 
%   stimulus and the EEG response.
% - In the future we will update this to accommodate both forward and 
%   backward time shifts (for instance, if the EEG anticipates stimulus
%   feature fluctuations).
filterLength=fsStim; 

%%% Now create the convolution matrix
% - This is the matrix with all the shifted versions of the stimulus
%   feature.
% OLD: sampleFeatureConvolution=tplitz(sampleFeature,filterLength);

% Inputs to createFeatureToeplitzMatrix:
% - Input 1 is the stimulus feature vector
% - Input 2 is the number of delayed versions we want to add. For now we
%   want fsStim instances total, included the non-delayed version, so we
%   specify filterLength-1 delays.
% - Input 3 is whether to add the intercept column at the end
zStimRMSMatrix = createFeatureToeplitzMatrix(zStimRMS, filterLength-1, 1);
zStimFluxMatrix = createFeatureToeplitzMatrix(zStimFlux, filterLength-1, 1);
zStimPC1Matrix = createFeatureToeplitzMatrix(zStimPC1, filterLength-1, 1);

%% Specify regularization parameters
% We are ready to call the core (CCA) function, but first we need to set
% some regularization parameters, which tell the CCA how many dimensions we
% should keep in both the stimulus and the EEG data
% - Blair note: Maybe leave these regularization variables as is for now, 
%   and we can see whether changing them makes a difference in the future 
%   (with RCA we haven't tried to tune these much)
% - Smaller number means stronger regularization

% How strongly to regularize the stimulus
Kx=7; % Starting value = 7

% How strongly to regularize the EEG
Ky=7; % Starting value = 7

%% Call the core CCA function

%%% "X" refers to stimulus things, and "Y" refers to response things

%%% As a reminder, here are the current possible inputs
% - eegDC: Downsampled, DC-corrected [space x time x trial] EEG matrix
% - zStimRMSMatrix: Convolution matrix of RMS feature
% - zStimFluxMatrix: Convolution matrix of Flux feature
% - zStimPC1Matrix: Convolution matrix of PC1 feature

%%% Reshape the EEG into a single space x concat time matrix
eegConcat = eegDC(:, :)'; % [concat time x space]

%%% Call the core function which computes the spatial and temporal filters
% - H: Temporal filters
% - W: Spatial filters
% - U: Temporally filtered stim feature
% - V: Spatially filtered EEG
% - A: Fwd model of spatial filters
% - dC: Eigenvalues, R: Struct of covariance matrices
% - R: Struct containing all the covariance matrices

%%% Specify which feature to input and repeat it nTrials times
featureUse = zStimPC1Matrix; 
fStr = 'PC1';
featureUseConcat = repmat(featureUse, nTrial, 1); % [concat time x lag]

[H, W, U, V, A, dC, R] = computeCCA(featureUseConcat, eegConcat, Kx, Ky); 

%% Compute SRC - multi trial

%%% For multiple trials, we could technically do one big correlation on 
% the per-component concatenated trials and stim features (the next code
% block). But we can also do the following (this code block):
% - Reshape the stim feature and EEG back to 3d [time x space x trial]
%   matrices
% - For each component, compute SRC for each trial
% - Report the average SRC across all trials

%%% Reshape filtered stim feature matrix
% Note that each copy (trial) of the stim feature is the same since we were
% multiplying a repeated single input feature matrix by the same temporal
% filter
nComp = min(size(U, 2), size(V, 2)); % Number of components to work with
U3d = reshape(U', [nComp, nTime, nTrial]); % [space x time x trial[
U3d = permute(U3d, [2 1 3]); % [time x space x trial]

%%% Reshape filtered EEG matrix
V3d = reshape(V', [nComp, nTime, nTrial]); % [space x time x trial]
V3d = permute(V3d, [2 1 3]); % [time x space x trial]

%%% Now iterate through the trials AND the components
% - Looping for clarity; can probably vectorize for speed later

allRho = nan(nTrial, nComp); % We will store all rho values in here
for i = 1:nTrial

    % Get the current U (filtered stim feature) and V (filtered EEG)
    thisU = U3d(:, :, i); % Time x component
    thisV = V3d(:, :, i); % Time x component

    for j = 1:nComp
        
        allRho(i, j) = corr(thisU(:, j), thisV(:, j));

    end

    clear this*

end

% Now take the mean across the trial dimension to get the mean SRC for each
% component
mean(allRho, 1)'

%% Compute SRC - single trial

%%% NEW 2/3/2026: This is now separate from the core CCA function
% This implementation keeps the outputs of the CCA function as one big
% filtered feature and one big filtered stimulus. 

% For how many components can we calculate SRC? 
nComp = min(size(U, 2), size(V, 2));
rhos = nan(nComp, 1);

% Do the correlation
for i = 1:nComp
    rhos(i) = corr(U(:,i), V(:,i));
end

rhos

%% Visualize results

% TODO: Add option for consistent clim and ylim across plots

close all
nPlot = 3; 

figure(1)
for c=1:nPlot
    subplot(2,nPlot,c);
    
    %%% BK: Topomap of current CC (spatial filter of EEG)
    plotOnEgi_multi(A(:,c), 4);
    title(['Spatial response: comp ' num2str(c)]);
    colormap(jmaColors('coolhotcortex'))
    
    %%% BK: Temporal filter of current CC (of stim feature)
    subplot(2,nPlot,c+nPlot);
    plot(H(1:end-1,c),'k');
    title(['Temporal response: comp. ' num2str(c)]);
end
sgtitle(fStr)

return

%% Visualize results - not working right now with EGI data

% 2/24/26: I'm having trouble getting the general-purpose topoplot function
% working with EGI sfp files right now - please use preceding block and
% plotOnEgi_multi() function for right now (first input is values to be
% visualized, second input is marker size on the head map)

%%% BK: Unlike in the core RCA function, the data are not visualized in the
%%% function that does the calculations but rather separately as shown
%%% below

% now let's look at the scalp maps of the first three EEG components that
% the stim2eeg found

%%% NEW 2/3/2026: When visualizing the temporal filters, omit the last
%%% value since that is the intercept

%%% TODO: 
% - Visualize topos on shared cLim and temporal filters on shared yLim
% - Change temporal filter xaxis to time (msec or sec) rather than samp

close all
nPlot = 3; 

%%% Figure 1: Close to original demo.m figure, except omitting last sample
%%% of the temporal filters (i.e., the intercept values)
figure(1);
locfile='GSN-HydroCel-125.sfp'; % EEGLAB-style location file for rendering scalp maps
for c=1:nPlot
    subplot(2,nPlot,c);
    
    %%% BK: Topomap of current CC (spatial filter of EEG)
    topoplot_new(A(:,c),locfile,'electrodes','off','numcontour',0,'plotrad',0.7);
    title(['Spatial response: comp ' num2str(c)]);
    colormap('jet')
    
    %%% BK: Temporal filter of current CC (of stim feature)
    subplot(2,nPlot,c+nPlot);
    plot(H(1:end-1,c),'k');
    title(['Temporal response: comp. ' num2str(c)]);
end
sgtitle('Original plot (matches demo.m)')

%%% Figure 2: Plot topos on shared CLim, and line plots on shared yLim
APlot = A(:, 1:nPlot);
HPlot = H(1:end-1, 1:nPlot); % Also truncating last value (intercept)
AMax = max(abs(APlot(:))); % Get global abs max - will do symmetric CLim
HMax = max(abs(HPlot(:))); % Again global abs max for symmetric YLim
t = 1/fsStim * (0:size(HPlot,1)-1);
figure(2);
for c=1:nPlot
    subplot(2,nPlot,c);
    
    %%% BK: Topomap of current CC (spatial filter of EEG)
    topoplot_new(APlot(:,c),locfile,'electrodes','off','numcontour',0,'plotrad',0.7);
    title(['Spatial response: CC' num2str(c)]);
    colormap('jet')
    set(gca, 'CLim', [-AMax AMax])
    if c == 1
        tSize = get(gca, 'Position');
        h = colorbar; h.Location = 'westoutside';
        hTitle = get(h, 'Title'); 
        set(hTitle, 'String', 'A.U.')
        set(gca, 'Position', tSize)
    end
    
    %%% BK: Temporal filter of current CC (of stim feature)
    subplot(2,nPlot,c+nPlot);
    plot(t, HPlot(:,c),'k'); grid on
    xlabel('Time (sec)'); ylabel('Weight')
    title(['Temporal response: CC' num2str(c)]);
    set(gca, 'YLim', [-HMax HMax])
end
sgtitle('Updated plot (shared CLim and YLim)')



