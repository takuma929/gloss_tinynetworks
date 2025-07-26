%--------------------------------------------------------------------------
% Online Experimental Data Processing Script
%
% This script processes raw online participant data for the gloss judgement
% experiment. It loads response files for all observers, 
% applies two exclusion criteria (failed attention/catch trial,
% and abnormally short response times), computes derived gloss parameters
% (Pellacini's c and Ward's specularity), and links each response to its
% corresponding ground-truth metadata. The script outputs cleaned data
% structures (including only valid participants) for downstream analysis and
% saves them as 'onlineData.mat'.
%
% Key Steps:
%   - Load and parse participant response files
%   - Apply exclusion criteria:
%       (1) Failure of catch trial (attention check)
%       (2) Median response time less than 3 s
%   - Convert raw slider responses to perceptual gloss parameters
%   - Organize all responses and ground-truth labels for each valid observer
%   - Save the cleaned dataset for further analysis
%
% Outputs:
%   - 'data/onlineData.mat': Cleaned online experiment dataset
%
% Author: TM, 2025
%--------------------------------------------------------------------------

disp('Processing online experiment data...')
clearvars; close all; % Clean workspace and close all figures

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

%% ------------------- Load Experiment Parameters and Files ---------------------

% Load table summarizing ground-truth for all stimuli and conditions
table_gt = readtable(fullfile('data', 'onlineExp_condition_summary'),'VariableNamingRule', 'modify');

% Get list of all data files starting with 'ID'
files = dir(fullfile('onlineData', 'ID*'));
nFiles = length(files);

% Set up parameters for gloss computation
resolution = 52;
L = 50;
cmin = 0; cmax = 0.1487; % Pellacini's c bounds
c = linspace(cmin, cmax, resolution);       % Pellacini's c
ps = (c + (0.5^3 / 2)^(1/3)).^3 - 0.5^3 / 2; % Ward's specularity

% Trial indices for two repetitions
trialN_main(:, 1) = 9:93;    % Repetition 1
trialN_main(:, 2) = 94:178;  % Repetition 2

%% ------------------- Initialize Variables ---------------------

cnt = 1;                       % Counter for valid participants
cnt_passedObserver = 1;        % Counter for passed observers
data = struct();               % Store participant data
gt = struct();                 % Store ground-truth per observer
passdObsList = [];             % Indices of passed observers
medianresponseTime = {};       % Median RT storage
failedObsFileList = {};        % Filenames of failed observers

%% ------------------- Process Each Observer ---------------------

for N = 1:nFiles
    fname = files(N).name;
    table = readtable(fullfile('onlineData', fname),'VariableNamingRule', 'modify');
    
    % Store basic metadata
    data(N).age        = table.age(1);
    data(N).gender     = table.biologicalSex(1);
    data(N).date       = table.date(1);
    data(N).OS         = table.OS(1);
    data(N).frameRate  = table.frameRate(1);
    data(N).groupN     = table.condition(8);

    % Trial-by-trial data for both repetitions
    for repetitionN = 1:2
        trialIdx = trialN_main(:, repetitionN);
        imgN_temp = table.imgN(trialIdx);
        [~, sortedId] = sort(imgN_temp);

        trial_valid = sortedId(2:end);   % Ignore first as catch trial
        trial_catch = sortedId(1);       % First is catch trial

        data(N).imgN(:, repetitionN)           = table.imgN(trialIdx(trial_valid));
        data(N).responseTime(:, repetitionN)   = table.key_resp_rt(trialIdx(trial_valid));
        data(N).response_rawN(:, repetitionN)  = table.slider_response(trialIdx(trial_valid));
        data(N).catch_response(:, repetitionN) = table.slider_response(trialIdx(trial_catch));
        data(N).start_specularN(:, repetitionN)= table.start_specularN(trialIdx(trial_valid));

        % Replace missing responses with start_specularN
        missing = isnan(data(N).response_rawN(:, repetitionN));
        data(N).response_rawN(missing, repetitionN) = data(N).start_specularN(missing, repetitionN);
    end

    % Convert response index to Pellacini's c and Ward's specularity
    data(N).response_Pellacini_c = c(data(N).response_rawN);
    data(N).response_specularity = ps(data(N).response_rawN);

    % ---- Retrieve corresponding ground-truth values for each image ----
    for imgN = 1:72
        idx = (table_gt.groupN == data(N).groupN) & (table_gt.testimgN == imgN);
        gt(N).lightprobeN(imgN, 1) = table_gt.lightprobeN(idx);
        gt(N).objN(imgN, 1)        = table_gt.objN(idx);
        gt(N).viewPointN(imgN, 1)  = table_gt.viewPointN(idx);
        gt(N).specularity(imgN, 1) = table_gt.specularity(idx);
        gt(N).Pellacini_c(imgN, 1) = table_gt.Pellacini_c(idx);
        gt(N).hue(imgN, 1)         = table_gt.hue(idx);
        gt(N).chroma(imgN, 1)      = table_gt.chroma(idx);
    end

    % ---- Load ground-truth for 12 common images from previous experiment ----
    load(fullfile('data', 'groundtruth_prev_exp'), 'GroundTruth');
    cnt_img = 0;
    for imgN = 3:3:38
        cnt_img = cnt_img + 1;
        imgIdx = 72 + cnt_img;
        gt(N).lightprobeN(imgIdx, 1) = GroundTruth.lightProbeN(imgN);
        gt(N).objN(imgIdx, 1)        = nan;
        gt(N).viewPointN(imgIdx, 1)  = nan;
        gt(N).specularity(imgIdx, 1) = GroundTruth.specularity(imgN);
        gt(N).Pellacini_c(imgIdx, 1) = GroundTruth.Pellacini_c(imgN);
        gt(N).hue(imgIdx, 1)         = GroundTruth.hue(imgN);
        gt(N).chroma(imgIdx, 1)      = GroundTruth.chroma(imgN);
    end

    % ---- Observer Exclusion Criteria ----
    % Catch trial (should sum to 2 for both repetitions)
    sum_catch = sum(data(N).catch_response(:));
    if sum_catch <= 4
        data(N).attentionPass = 1;
    else
        data(N).attentionPass = 0;
        failedObsFileList{end+1,1} = fname; % Store failed filename
    end

    % Median response time criterion
    if data(N).attentionPass
        medianRT = median(data(N).responseTime(:));
        medianresponseTime{cnt, 1} = medianRT;
        if medianRT > 3  % Criterion based on reference RT
            passdObsList(cnt_passedObserver,1) = N;
            cnt_passedObserver = cnt_passedObserver + 1;
        else
            failedObsFileList{end+1,1} = fname; % Also failed due to RT
        end
        cnt = cnt + 1;
    end
end

% Only passed observers are included in final datasets
data = data(passdObsList); 
gt = gt(passdObsList);

% Optionally save cleaned data
fname = fullfile('data','onlineData.mat');
save(fname, 'data', 'gt');
disp(['Saved: ',fname])

