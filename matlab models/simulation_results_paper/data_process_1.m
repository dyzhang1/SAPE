% BLER per Subband under Deep Fading - Predistortion vs Non-Predistortion
% This script compares BLER per subband for predistortion vs non-predistortion
% configurations at a specified SNR and CQI under deep fading conditions.

%% 1. Load simulation results for predistortion and non-predistortion
pred_data = load('predistortion_results.mat');
nonpred_data = load('nonpredistortion_results.mat');
% (Adjust file names and paths as needed. The .mat files are assumed to contain 
% variables such as SNR_values, CQI_values, and subband_BLER arrays.)

%% 2. Define parameters: target SNR and CQI for analysis
targetSNR_dB = 4;   % SNR value in dB (e.g., 2 dB)
targetCQI    = 9;   % CQI level (e.g., 3 corresponds to CQI3)

%% 3. Extract BLER values per subband at the specified SNR and CQI
% Find index of the SNR closest to targetSNR_dB in the loaded data
[~, snrIndex] = min(abs(pred_data.SNR_values - targetSNR_dB));
% If data contains multiple CQI values, find the index for targetCQI
if isfield(pred_data, 'CQI_values')
    cqiIndex = find(pred_data.CQI_values == targetCQI);
    bler_pred_subbands = pred_data.subband_BLER(:, cqiIndex, snrIndex);
    bler_non_subbands  = nonpred_data.subband_BLER(:, cqiIndex, snrIndex);
else
    % If each file is for a single CQI, subband_BLER may be 2D (subband x SNR)
    bler_pred_subbands = pred_data.subband_BLER(:, snrIndex);
    bler_non_subbands  = nonpred_data.subband_BLER(:, snrIndex);
end
% Ensure the BLER data is in column-vector form for plotting
bler_pred_subbands = bler_pred_subbands(:);
bler_non_subbands  = bler_non_subbands(:);

%% 4. Plot a grouped bar chart for BLER per subband
figure;
% Combine data for both scenarios into a matrix for grouped bar plot
bler_matrix = [bler_non_subbands, bler_pred_subbands];
bar(bler_matrix, 'grouped');
% Using a matrix with two columns creates side-by-side bars for each subband.

% Label the axes
xlabel('Subband Index');
ylabel('Block Error Rate (BLER)');
% Optionally, add a title for clarity
title(sprintf('Subband BLER at %g dB (CQI %d)', targetSNR_dB, targetCQI));

% Customize x-axis ticks to label each subband
numSubbands = length(bler_pred_subbands);
xticks(1:numSubbands);
xticklabels(1:numSubbands);

% Add a legend to distinguish the scenarios
leg = legend('No Predistortion', 'With Predistortion', 'Location', 'best');
set(leg, 'FontSize', 14);  % increase legend font size
% (Optional: set(leg, 'Box', 'off') to remove legend border)

%% 5. Format the plot for publication-quality visuals
set(gca, 'FontSize', 14);                  % increase tick label font size
xlabel('Subband Index', 'FontSize', 16);   % increase x-axis label font size
ylabel('Block Error Rate (BLER)', 'FontSize', 16);  % increase y-axis label font size
% (You can also adjust font name, line width, etc., as needed.)
