% This code creates a full virtual cohort for the model
% We generate 16 here concatanete them with other 3 (r1,r2,K1) coming from 
% "generate_vp_control.m" as discussed therein. We use the following bounds
% for all 16 parameters
                % test1--[0.25*base,1.50*base]
                % test2--[0.70*base,1.30*base]

clc;clear all;
% Define baseline values
baseline_l = 1.56;
baseline_dC = 0.41;
baseline_muC = 0.6;
baseline_s = 3.05e-1;
baseline_gammaC = 2.93e-2;
baseline_k = 2.019e-7;
baseline_wC = 3e-5;
baseline_K2 = 1.65e3;
baseline_dB = 0.3;
baseline_muB = 0.89;
baseline_b = 1.4e-3;
baseline_gammaB = 7e-3;
baseline_wB = 3.42e-6;


% Load vp_control data, this data is fixed, we populate it with additional
% 16 parameters.
vp_control = table2array(readtable("meta_data/vp_control.csv"));
N = size(vp_control,1);

% define a convenient structure for for loop.
tests = struct('name', {'test1', 'test2'}, ...
               'low_factor', {0.25, 0.70}, ...
               'up_factor', {1.50, 1.30});

% Loop through each test case
for i = 1:length(tests)
    test = tests(i);

    % save in the relevant path, matlab's snytax is terrible anyway...
    vp_full = fullfile('meta_data', test.name, ['vp_full_' test.name '.csv']);
    disp(vp_full)



    % delete it if exist, sometimes matlab persists
    if exist(vp_full, 'file') == 2
        delete(vp_full);
    end

    % Generate random values within the specified range for each variable
    l = random_between_range(test.low_factor * baseline_l, test.up_factor * baseline_l, N);
    dC = random_between_range(test.low_factor * baseline_dC, test.up_factor * baseline_dC, N);
    muC = random_between_range(test.low_factor * baseline_muC, test.up_factor * baseline_muC, N);
    s = random_between_range(test.low_factor * baseline_s, test.up_factor * baseline_s, N);
    gammaC = random_between_range(test.low_factor * baseline_gammaC, test.up_factor * baseline_gammaC, N);
    k = random_between_range(test.low_factor * baseline_k, test.up_factor * baseline_k, N);
    wC = random_between_range(test.low_factor * baseline_wC, test.up_factor * baseline_wC, N);
    K2 = random_between_range(test.low_factor * baseline_K2, test.up_factor * baseline_K2, N);
    dB = random_between_range(test.low_factor * baseline_dB, test.up_factor * baseline_dB, N);
    muB = random_between_range(test.low_factor * baseline_muB, test.up_factor * baseline_muB, N);
    b = random_between_range(test.low_factor * baseline_b, test.up_factor * baseline_b, N);
    gammaB = random_between_range(test.low_factor * baseline_gammaB, test.up_factor * baseline_gammaB, N);
    wB = random_between_range(test.low_factor * baseline_wB, test.up_factor * baseline_wB, N);

    % Concatenate the data
    full_data = [vp_control, l, dC, muC, s, gammaC, k, wC, K2, dB, muB, b, gammaB, wB];
    headers = {'r1', 'r2', 'K1', 'l', 'dC', 'muC', 's', 'gammaC', 'k', 'wC', 'K2', ...
               'dB', 'muB', 'b', 'gammaB', 'wB'};

    % Convert to table and write to CSV
    T = array2table(full_data, 'VariableNames', headers);
    writetable(T, vp_full);
end

% Function to generate random values within a specified range
function z = random_between_range(lower_bound, upper_bound, N)
    z = lower_bound + (upper_bound - lower_bound) * rand(1, N)';
end













