% Modified main.m for Monthly to Daily Disaggregation
% Purpose: Convert monthly reports to daily estimates while maintaining the experimental structure
% Period covered: January 2020 - November 2024

addpath('HFusion')
addpath('AresFilter')
addpath('data')
clc; clear;

%% 1. Define monthly reports (Jan 2020 - Nov 2024)
reports = [
    1, 31, 43787; % Jan 2020
    32, 60, 40427; % Feb 2020 (29 days, leap year)
    61, 91, 39398; % Mar 2020
    92, 121, 25369; % Apr 2020
    122, 152, 19049; % May 2020
    153, 182, 26593; % Jun 2020
    183, 213, 26878; % Jul 2020
    214, 244, 23179; % Aug 2020
    245, 274, 24317; % Sep 2020
    275, 305, 20644; % Oct 2020
    306, 335, 21124; % Nov 2020
    336, 366, 16747; % Dec 2020
    367, 397, 35053; % Jan 2021
    398, 425, 31532; % Feb 2021 (28 days)
    426, 456, 39848; % Mar 2021
    457, 486, 38948; % Apr 2021
    487, 517, 34749; % May 2021
    518, 547, 37652; % Jun 2021
    548, 578, 25075; % Jul 2021
    579, 609, 30670; % Aug 2021
    610, 639, 39631; % Sep 2021
    640, 670, 42452; % Oct 2021
    671, 700, 45619; % Nov 2021
    701, 731, 44017; % Dec 2021
    732, 762, 66025; % Jan 2022
    763, 790, 48431; % Feb 2022 (28 days)
    791, 821, 56857; % Mar 2022
    822, 851, 52675; % Apr 2022
    852, 882, 54903; % May 2022
    883, 912, 66734; % Jun 2022
    913, 943, 63507; % Jul 2022
    944, 974, 67759; % Aug 2022
    975, 1004, 69874; % Sep 2022
    1005, 1035, 69597; % Oct 2022
    1036, 1065, 68806; % Nov 2022
    1066, 1096, 56626; % Dec 2022
    1097, 1127, 83897; % Jan 2023
    1128, 1155, 73983; % Feb 2023 (28 days)
    1156, 1186, 77314; % Mar 2023
    1187, 1216, 56881; % Apr 2023
    1217, 1247, 80216; % May 2023
    1248, 1277, 69611; % Jun 2023
    1278, 1308, 78302; % Jul 2023
    1309, 1339, 85188; % Aug 2023
    1340, 1369, 78406; % Sep 2023
    1370, 1400, 83093; % Oct 2023
    1401, 1430, 77180; % Nov 2023
    1431, 1461, 63072; % Dec 2023
    1462, 1492, 74898; % Jan 2024
    1493, 1521, 65049; % Feb 2024 (29 days, leap year)
    1522, 1552, 67988; % Mar 2024
    1553, 1582, 66322; % Apr 2024
    1583, 1613, 75951; % May 2024
    1614, 1643, 65180; % Jun 2024
    1644, 1674, 74885; % Jul 2024
    1675, 1705, 71239; % Aug 2024
    1706, 1735, 70057; % Sep 2024
    1736, 1766, 77470; % Oct 2024
    1767, 1796, 70135; % Nov 2024
];

% Get total number of days
num_days = reports(end, 2);

%% 2. Parameter Optimization
fprintf('Starting parameter optimization...\n');
[best_params, best_score] = parameter_optimization(reports, num_days);

% Use the best parameters found
lambdas = best_params.lambda;
alpha = best_params.alpha;
gama = best_params.gama;
iterationTime = best_params.iterationTime;
stop_ratios = best_params.stop_ratio;

fprintf('\nUsing optimized parameters:\n');
fprintf('Lambda: %.1f\n', lambdas);
fprintf('Alpha: %.1f\n', alpha);
fprintf('Gama: %.1f\n', gama);
fprintf('Iteration Time: %d\n', iterationTime);
fprintf('Stop Ratio: %d\n', stop_ratios);

%% 3. First Phase: Initial reconstruction using H-Fusion
fprintf('Starting H-Fusion reconstruction...\n');

% Get constraint matrix from the reports
[A, y] = rep_constraint_equations_full(reports, num_days);

% Save the constraint equations
all_results.constraints = struct('A', A, 'y', y);

% Create initial reconstruction using original reports
Out_original = struct();
Out_original(1).muvars = [0, 0];  % Special indicator for original reports
Out_original(1).A = A;
Out_original(1).y = y;

% Create a temporary events vector just for sp_reconstruct
temp_events = zeros(num_days, 1);
[recon_events, recon_error, reconstruction_param, M] = sp_reconstruct(A, y, lambdas, temp_events, alpha);
Out_original(1).x_reconstr = recon_events(:, 1, 1);  % Use first lambda and alpha values
Out_original(1).x_error = recon_error;
Out_original(1).Matrix = M;
Out_original(1).sp_params = reconstruction_param;
[Out_original(1).error, Out_original(1).minIdx] = min(recon_error(:));

% Save H-Fusion initial results
all_results.hfusion_initial = struct('recon_events', recon_events, 'recon_error', recon_error, ...
                                    'reconstruction_param', reconstruction_param, 'M', M, ...
                                    'Out_original', Out_original);

% Combine original reports results
Out = Out_original;
all_results.hfusion_combined = Out;

%% 4. Second Phase: Enhance reconstruction using Annihilating Filter
fprintf('Starting ARES filter processing...\n');

% Try different annihilating filter ratios
fprintf('Testing with stop ratio = %.2f\n', stop_ratios(1));

% Initialize error metrics arrays
error_metrics_non_int = struct('monthly_avg', zeros(1,iterationTime), ...
                             'neighboring', zeros(1,iterationTime), ...
                             'weekly', zeros(1,iterationTime), ...
                             'monthly', zeros(1,iterationTime));
error_metrics_int = struct('monthly_avg', zeros(1,iterationTime), ...
                         'neighboring', zeros(1,iterationTime), ...
                         'weekly', zeros(1,iterationTime), ...
                         'monthly', zeros(1,iterationTime));

% Save ARES parameters
all_results.ares_params = struct('stop_ratios', stop_ratios);

% For each ratio, run the iteration process
all_results.ares_results = struct();
for i = 1:length(stop_ratios)
    fprintf('\n--- Testing Stop Ratio %d of %d (value: %d) ---\n', i, length(stop_ratios), stop_ratios(i));
    
    % Run iterative refinement with annihilating filter
    fprintf('Running %d iterations...\n', iterationTime);
    [Inc_A, Out_A1, Out_AL] = iteration(Out, num_days, iterationTime, gama, stop_ratios(i));
    
    % Save results for this ratio
    fprintf('Saving results for ratio %d...\n', stop_ratios(i));
    ratio_results = struct('Inc_A', Inc_A, 'Out_A1', Out_A1, 'Out_AL', Out_AL);
    all_results.ares_results(i).ratio = stop_ratios(i);
    all_results.ares_results(i).data = ratio_results;
    
    % Process results (since we only have one ratio)
    fprintf('\nProcessing results...\n');
    daily_estimates = Out_AL(1).x_reconstr;
    
    % Get H-Fusion results for comparison
    hfusion_daily = Out_original(1).x_reconstr;
    
    % Calculate smoothness metrics for both methods
    fprintf('Calculating smoothness metrics...\n');
    hfusion_smoothness = calculate_smoothness(hfusion_daily);
    ares_smoothness = calculate_smoothness(daily_estimates);
    
    % Calculate daily variations
    hfusion_variations = calculate_daily_variations(hfusion_daily);
    ares_variations = calculate_daily_variations(daily_estimates);
    
    % Store comparison metrics
    all_results.method_comparison = struct(...
        'hfusion', struct(...
            'smoothness', hfusion_smoothness, ...
            'daily_variations', hfusion_variations), ...
        'ares', struct(...
            'smoothness', ares_smoothness, ...
            'daily_variations', ares_variations));
    
    % Plot comparison
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot daily estimates
    subplot(2,1,1);
    plot(1:num_days, hfusion_daily, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(1:num_days, daily_estimates, 'r-', 'LineWidth', 1.5);
    title('Daily Estimates Comparison');
    xlabel('Day');
    ylabel('Number of Cases');
    legend('H-Fusion', 'ARES');
    grid on;
    
    % Plot monthly averages
    subplot(2,1,2);
    monthly_hfusion = zeros(size(reports,1), 1);
    monthly_ares = zeros(size(reports,1), 1);
    for j = 1:size(reports,1)
        monthly_hfusion(j) = mean(hfusion_daily(reports(j,1):reports(j,2)));
        monthly_ares(j) = mean(daily_estimates(reports(j,1):reports(j,2)));
    end
    bar([monthly_hfusion, monthly_ares]);
    title('Monthly Average Comparison');
    xlabel('Month');
    ylabel('Average Daily Cases');
    legend('H-Fusion', 'ARES');
    grid on;
    
    % Save comparison plot
    saveas(gcf, 'method_comparison.png');
    
    % Apply integer constraints
    fprintf('Applying integer constraints...\n');
    % First get non-integer version that matches monthly totals exactly
    non_integer_daily = enforce_integer_constraints(daily_estimates, reports, false);
    % Then get integer version that also matches monthly totals
    integer_daily = enforce_integer_constraints(daily_estimates, reports, true);
    
    % Calculate error metrics for both versions
    fprintf('Calculating error metrics...\n');
    error_metrics_non_int = calculate_error_metrics(non_integer_daily, reports);
    error_metrics_int = calculate_error_metrics(integer_daily, reports);
    
    % Store both versions and their error metrics
    all_results.final_daily_estimates = struct('non_integer', non_integer_daily, ...
                                             'integer', integer_daily, ...
                                             'error_metrics_non_int', error_metrics_non_int, ...
                                             'error_metrics_int', error_metrics_int);
    
    % Calculate monthly sums to verify constraints for both versions
    fprintf('Calculating monthly sums for verification...\n');
    monthly_sums_non_int = zeros(size(reports,1), 1);
    monthly_sums_int = zeros(size(reports,1), 1);
    
    for j = 1:size(reports,1)
        monthly_sums_non_int(j) = sum(non_integer_daily(reports(j,1):reports(j,2)));
        monthly_sums_int(j) = sum(integer_daily(reports(j,1):reports(j,2)));
    end
    
    % Calculate constraint validation metrics for both versions
    fprintf('Computing constraint validation metrics...\n');
    constraint_error_non_int = reports(:,3) - monthly_sums_non_int;
    constraint_error_int = reports(:,3) - monthly_sums_int;
    
    relative_error_non_int = abs(constraint_error_non_int) ./ reports(:,3) * 100;
    relative_error_int = abs(constraint_error_int) ./ reports(:,3) * 100;
    
    all_results.constraint_validation = struct(...
        'non_integer', struct(...
            'absolute_error', constraint_error_non_int, ...
            'relative_error', relative_error_non_int, ...
            'mean_abs_error', mean(abs(constraint_error_non_int)), ...
            'mean_rel_error', mean(relative_error_non_int)), ...
        'integer', struct(...
            'absolute_error', constraint_error_int, ...
            'relative_error', relative_error_int, ...
            'mean_abs_error', mean(abs(constraint_error_int)), ...
            'mean_rel_error', mean(relative_error_int)));
    
    fprintf('✓ Validation metrics computed\n');
    
    % Display error metrics
    fprintf('\nError Metrics Summary:\n');
    fprintf('Non-Integer Version:\n');
    fprintf('- Monthly Average ERROR: %.4f\n', error_metrics_non_int.monthly_avg);
    fprintf('- Neighboring Days ERROR: %.4f\n', error_metrics_non_int.neighboring);
    fprintf('- Weekly Pattern ERROR: %.4f\n', error_metrics_non_int.weekly);
    fprintf('- Monthly Pattern ERROR: %.4f\n', error_metrics_non_int.monthly);
    
    fprintf('\nInteger Version:\n');
    fprintf('- Monthly Average ERROR: %.4f\n', error_metrics_int.monthly_avg);
    fprintf('- Neighboring Days ERROR: %.4f\n', error_metrics_int.neighboring);
    fprintf('- Weekly Pattern ERROR: %.4f\n', error_metrics_int.weekly);
    fprintf('- Monthly Pattern ERROR: %.4f\n', error_metrics_int.monthly);
    
    % Display method comparison metrics
    fprintf('\nMethod Comparison Summary:\n');
    fprintf('H-Fusion:\n');
    fprintf('- Smoothness Score: %.4f\n', hfusion_smoothness);
    fprintf('- Average Daily Variation: %.4f\n', mean(hfusion_variations));
    
    fprintf('\nARES:\n');
    fprintf('- Smoothness Score: %.4f\n', ares_smoothness);
    fprintf('- Average Daily Variation: %.4f\n', mean(ares_variations));
end

% Plot error metrics comparison
figure;
error_types = {'Monthly Avg', 'Neighboring', 'Weekly', 'Monthly'};
non_int_errors = [error_metrics_non_int.monthly_avg, error_metrics_non_int.neighboring, ...
                 error_metrics_non_int.weekly, error_metrics_non_int.monthly];
int_errors = [error_metrics_int.monthly_avg, error_metrics_int.neighboring, ...
             error_metrics_int.weekly, error_metrics_int.monthly];

bar([non_int_errors; int_errors]');
title('Error Metrics Comparison');
xlabel('Error Type');
ylabel('Error Value');
set(gca, 'XTickLabel', error_types);
legend('Non-Integer', 'Integer');
grid on;
saveas(gcf, 'error_metrics_comparison.png');

%% Save all variables to .mat file
save('monthly_to_daily_results.mat', 'all_results');
save('complete_workspace.mat');

% Create a detailed summary for easy reference
fid = fopen('summary_report.txt', 'w');
fprintf(fid, 'Monthly to Daily Disaggregation Summary\n');
fprintf(fid, '=====================================\n\n');
fprintf(fid, 'Parameters Used:\n');
fprintf(fid, '- Lambda values: '); fprintf(fid, '%.2f ', lambdas); fprintf(fid, '\n');
fprintf(fid, '- Alpha: %.2f\n', alpha);
fprintf(fid, '- Gamma: %.2f\n', gama);
fprintf(fid, '- Iteration Time: %d\n', iterationTime);
fprintf(fid, '- Stop Ratios: '); fprintf(fid, '%d ', stop_ratios); fprintf(fid, '\n\n');

fprintf(fid, 'Method Comparison:\n');
fprintf(fid, '-----------------\n');
fprintf(fid, 'H-Fusion:\n');
fprintf(fid, '- Smoothness Score: %.4f\n', hfusion_smoothness);
fprintf(fid, '- Average Daily Variation: %.4f\n', mean(hfusion_variations));

fprintf(fid, '\nARES:\n');
fprintf(fid, '- Smoothness Score: %.4f\n', ares_smoothness);
fprintf(fid, '- Average Daily Variation: %.4f\n', mean(ares_variations));

fprintf(fid, '\nError Metrics Results:\n');
fprintf(fid, '--------------------\n');
fprintf(fid, 'Non-Integer Version:\n');
fprintf(fid, '- Monthly Average ERROR: %.4f\n', error_metrics_non_int.monthly_avg);
fprintf(fid, '- Neighboring Days ERROR: %.4f\n', error_metrics_non_int.neighboring);
fprintf(fid, '- Weekly Pattern ERROR: %.4f\n', error_metrics_non_int.weekly);
fprintf(fid, '- Monthly Pattern ERROR: %.4f\n', error_metrics_non_int.monthly);

fprintf(fid, '\nInteger Version:\n');
fprintf(fid, '- Monthly Average ERROR: %.4f\n', error_metrics_int.monthly_avg);
fprintf(fid, '- Neighboring Days ERROR: %.4f\n', error_metrics_int.neighboring);
fprintf(fid, '- Weekly Pattern ERROR: %.4f\n', error_metrics_int.weekly);
fprintf(fid, '- Monthly Pattern ERROR: %.4f\n', error_metrics_int.monthly);

fprintf(fid, '\nConstraint Validation:\n');
fprintf(fid, '--------------------\n');
fprintf(fid, 'Non-Integer:\n');
fprintf(fid, 'Mean Absolute Error: %.6f\n', all_results.constraint_validation.non_integer.mean_abs_error);
fprintf(fid, 'Mean Relative Error: %.6f%%\n', all_results.constraint_validation.non_integer.mean_rel_error);
fprintf(fid, 'Integer:\n');
fprintf(fid, 'Mean Absolute Error: %.6f\n', all_results.constraint_validation.integer.mean_abs_error);
fprintf(fid, 'Mean Relative Error: %.6f%%\n', all_results.constraint_validation.integer.mean_rel_error);

fclose(fid);

fprintf('Process complete. All results saved to .mat files.\n');
fprintf('To load the results later, use: load(''monthly_to_daily_results.mat'')\n');
fprintf('To load the entire workspace, use: load(''complete_workspace.mat'')\n');
fprintf('A summary report has been saved to: summary_report.txt\n');
