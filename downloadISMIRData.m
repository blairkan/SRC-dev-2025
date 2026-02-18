% downloadISMIRData.m
% ----------------------------------------
% Blair - 2/17/2026
%
% This script downloads ready-to-use EEG data from the Gang et al. (2017)
% ISMIR paper into the "Data" folder of the local instance of the repo if
% the files are not already there.

% History
% 2/17/26: Adapted from MatClassRSA, illustrative_0_downloadExampleData.m

clear all; close all; clc

%%% Get path(s) to the "Data" folder in the SRC-dev-2025 repo.
dataPath = what(['SRC-dev-2025' filesep 'Data']);

% We need exactly one path to check and possibly download data into.
if length(dataPath) == 0
    error('"Data" folder not found. Make sure that SRC-dev-2025 is on your local machine and that the entire repo is added to the Matlab path.')
elseif length(dataPath) > 1
    dataPath = dataPath(1);
    warning(['More than one instance of "SRC-dev-2025/Data" was found on the local machine. If data files need to be downloaded, they will be downloaded to the first indexed directory: ' newline dataPath.path])
end

%%% Create the struct array with URL and filename information

% Calls a helper function in this script
INFO = createInfoStruct;

%%% Specify which files to check for and download if needed

% Default vector (all files)
fileIdx = 1:length(INFO);
nFilesAll = length(fileIdx); % For printing messages later

%%%%%%%%%%%%%%%%%%%%% OPTIONAL USER SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%

% If the user does not need to check for/download all files in the list,
% they can specify a subset of file numbers here, which will overwrite
% the default list. Otherwise, comment out the line below.

% fileIdx = [1 2 3 4 5 6 7 8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check for each data file and download if not already there
nFilesDownloaded = 0;
nFiles = length(fileIdx);

for f = 1:nFiles
    currIdx = fileIdx(f);
    currURL = INFO(currIdx).url;
    currFN = INFO(currIdx).fn;
    disp(['File ' num2str(currIdx) ' of ' num2str(nFilesAll)])
    if ~exist(fullfile(dataPath.path, currFN))

        disp(['Downloading ' currFN ' to ' dataPath.path newline])
        websave(fullfile(dataPath.path, currFN), currURL);
        nFilesDownloaded = nFilesDownloaded + 1;

    else
        disp([currFN ' already found in ' dataPath.path newline])
    end

end

disp(['\ * \ * Run complete. ' num2str(nFilesDownloaded) ' file(s) downloaded. * / * /'])

clear

%% Helper function

function INFO = createInfoStruct()

% This is a struct array. Each element has the (1) download URL and
% (2) name of to-be-saved file for the items in the SDR repository.

% Index 1: Song 21 data (a and b combined)
INFO(1).url = 'https://www.dropbox.com/scl/fi/ihnk69qhm8dlqnajfport/song21_all_Imputed.mat?rlkey=sfej3vrec1wdm2y2wlyrehfly&st=kckv682n&dl=1';
INFO(1).fn = 'song21_all_Imputed.mat';

end