% demoComputeCCA.m
% ------------------
% Blair dev version 2026
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
% 2/3/2026 - adapted from demo.m
% thanks to Paul DeGuzman for providing data
% (c) Jacek P. Dmochowski, 2019
% jdmochowski@ccny.cuny.edu

clear all; close all; clc

%% Load input data
addpath('..');
load('./sampleData.mat','sampleEEG','sampleFeature','fsEEG','fsStim');

%% Downsample EEG to fs of stim feature
% first downsample the EEG to the sampling rate of the stimulus
%%% BK: We would generally downsample to the lower sampling rate. We can also
% resample to the higher one but there will be no added information
sampleEEGdown=resample(sampleEEG,fsStim,fsEEG); 

%% z-score the stim feature
% normalize the feature because correlation doesn't care about scale
sampleFeature=zscore(sampleFeature);

%% Adjust length of stim feature to match length of EEG
% it's common for the downsampled EEG and the stimulus to have slightly
% different lengths, so we fix this here
%%% BK: Currently the procedure is to trim or zero pad the stim feature
%%% input to match the EEG, and never adjust the EEG

% check if feature missing samples
%%% BK: If the EEG is longer than the stim feature, zero-pad the end of the
%%% stim feature so that the lengths match
if size(sampleEEGdown,1)>numel(sampleFeature) 
    nMissing=size(sampleEEGdown,1)-numel(sampleFeature);
    sampleFeature=cat(1,sampleFeature,zeros(nMissing,1));
end

% check if feature has too many samples
%%% BK: If the EEG is shorter than the stimulus feature, truncate the
%%% stimulus feature to match the length of the EEG
if size(sampleEEGdown,1)<numel(sampleFeature)  
    sampleFeature=sampleFeature(1:size(sampleEEGdown,1));
end

%% Create matrix of time-shifted versions of stim feature
% now we need to create a convolution matrix from the one-dimensional
% feature time series.  this allows us to temporally filter the stimulus
% time series using a matrix-vector product.  tplitz.m is a function that
% creates the convolution matrix, so we call it here on the stimulus

% before calling the function, we need to specify how long we want the
% filter to be.  here we set this to one second worth of stimulus.  so
% we're looking back 1 second in time.  you can think of this as the
% maximum delay between the stimulus and the EEG response.
%%% BK: In the future we will update this to accommodate both forward and
%%% backward time shifts
filterLength=fsStim; 

% now create the convolution matrix
%%% BK: this is the matrix with all the shifted versions of the stimulus
%%% feature. 
sampleFeatureConvolution=tplitz(sampleFeature,filterLength);

%% Specify regularization parameters
% we are ready to call the core (CCA) function, but first we need to set
% some regularization parameters, which tell the CCA how many dimensions we
% should keep in both the stimulus and the EEG data

%%% BK: Maybe leave these regularization variables as is for now, and we
%%% can see whether changing them makes a difference in the future (with
%%% RCA we haven't tried to tune these much)
% how strongly to regularize the stimulus (small number means strong
% regularization)
Kx=7; % Starting value = 7

% how strongly to regularize the EEG
Ky=7; % Starting value = 7

%% Call the core CCA function

%%% BK: "X" refers to stimulus things, and "Y" refers to response things

%%% call the core function which correlates the stimulus with the EEG
% H: Temporal filters, W: Spatial filters, A: Fwd model of spatial filters
% U: Temporally filtered stim feature, V: Spatially filtered EEG
% dC: Eigenvalues, R: Struct of covariance matrices
[H, W, U, V, A, dC, R] = computeCCA(sampleFeatureConvolution,sampleEEGdown,Kx,Ky); 

%% Compute SRC

%%% NEW 2/3/2026: This is now separate from the core CCA function

% For how many components can we calculate SRC? 
nComp = min(size(U, 2), size(V, 2));
rhos = nan(nComp, 1);

% Do the correlation
for i = 1:nComp
    rhos(i) = corr(U(:,i), V(:,i));
end

rhos

%% Visualize results

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

figure(1);
locfile='BioSemi32.loc'; % EEGLAB-style location file for rendering scalp maps
for c=1:3
    subplot(2,3,c);
    
    %%% BK: Topomap of current CC (spatial filter of EEG)
    topoplot_new(A(:,c),locfile,'electrodes','off','numcontour',0,'plotrad',0.7);
    title(['Spatial response: comp ' num2str(c)]);
    colormap('jet')
    
    %%% BK: Temporal filter of current CC (of stim feature)
    subplot(2,3,c+3);
    plot(H(1:end-1,c),'k');
    title(['Temporal response: comp. ' num2str(c)]);
end
