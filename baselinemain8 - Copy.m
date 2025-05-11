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

% Create placeholder events vector for the full time period
num_days = reports(end, 2);
events = zeros(num_days, 1);

% Initialize events from reports
for i = 1:size(reports, 1)
    start_day = reports(i, 1);
    end_day = reports(i, 2);
    total_events = reports(i, 3);
    % Distribute events evenly across the month
    daily_events = total_events / (end_day - start_day + 1);
    events(start_day:end_day) = daily_events;
end

%% 2. Setup parameters
% H-FUSION Parameters
lambdas = [0.1, 1, 5, 10, 20];  % Multiple lambda values to test
alpha = 0.5;                    % Balance between smoothness and periodicity
gama = 0.5;                     % Parameter for cost function
iterationTime = 4;              % Number of iterations

% Configure different experimental setups for report aggregation tests
% Well use sliding windows of reports with different durations and overlaps
% For your monthly data use case, this is primarily for testing the robustness of the method
config_rep_dur = [15, 30, 60];  % Report durations: half-month, month, 2-months
config_rep_over = [1, 7, 15];   % Report overlaps: daily, weekly, bi-weekly
xdim = length(config_rep_dur);
ydim = length(config_rep_over);

% Create a struct to store all important variables and results
all_results = struct();
all_results.parameters = struct('lambdas', lambdas, 'alpha', alpha, 'gama', gama, ...
                               'iterationTime', iterationTime, 'config_rep_dur', config_rep_dur, ...
                               'config_rep_over', config_rep_over);
all_results.reports = reports;

%% 3. First Phase: Initial reconstruction using H-Fusion
fprintf('Starting H-Fusion reconstruction...\n');

% Get constraint matrix from the reports
[A, y] = rep_constraint_equations_full(reports, events);

% Save the constraint equations
all_results.constraints = struct('A', A, 'y', y);

% Create initial reconstruction using original reports
Out_original = struct();
Out_original(1).muvars = [0, 0];  % Special indicator for original reports
Out_original(1).A = A;
Out_original(1).y = y;

[recon_events, recon_error, reconstruction_param, M] = sp_reconstruct(A, y, lambdas, events, alpha);
Out_original(1).x_reconstr = recon_events(:, 1, 1);  % Use first lambda and alpha values
Out_original(1).x_error = recon_error;
Out_original(1).Matrix = M;
Out_original(1).sp_params = reconstruction_param;
[Out_original(1).error, Out_original(1).minIdx] = min(recon_error(:));

% Save H-Fusion initial results
all_results.hfusion_initial = struct('recon_events', recon_events, 'recon_error', recon_error, ...
                                    'reconstruction_param', reconstruction_param, 'M', M, ...
                                    'Out_original', Out_original);

% Also create experimental synthetic reports of different durations to test robustness
fprintf('Running H-Fusion with different report configurations...\n');
[Out, Out_LSQ] = hfusion(events, lambdas, alpha, config_rep_dur, config_rep_over);

% Save all H-Fusion configuration results
all_results.hfusion_configs = struct('Out', Out, 'Out_LSQ', Out_LSQ);

% Combine original reports results with synthetic test results
Out = [Out_original, Out];
all_results.hfusion_combined = Out;

%% 4. Second Phase: Enhance reconstruction using Annihilating Filter
fprintf('Starting ARES filter processing...\n');

% Find optimal stop ratio using selection_A
fprintf('Finding optimal stop ratio...\n');
[AUC, Out_A_optimal] = selection_A(Out, events, Out);  % Using Out as both input and comparison
optimal_stop_ratio = AUC.final_stop;
fprintf('Optimal stop ratio found: %.2f\n', optimal_stop_ratio);

% Try different annihilating filter ratios
stop_ratios = [optimal_stop_ratio];  % Use the optimal ratio found
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
    [Inc_A, Out_A1, Out_AL] = iteration(Out, events, iterationTime, gama, stop_ratios(i));
    
    % Save results for this ratio
    fprintf('Saving results for ratio %d...\n', stop_ratios(i));
    ratio_results = struct('Inc_A', Inc_A, 'Out_A1', Out_A1, 'Out_AL', Out_AL);
    all_results.ares_results(i).ratio = stop_ratios(i);
    all_results.ares_results(i).data = ratio_results;
    
    % Process results (since we only have one ratio)
    fprintf('\nProcessing results...\n');
    daily_estimates = Out_AL(1).x_reconstr;
    
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
% First save the main results structure
save('monthly_to_daily_results.mat', 'all_results');

% Also save individual important variables for direct access
save('complete_workspace.mat');  % This saves ALL variables in the workspace

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

fprintf(fid, 'Error Metrics Results:\n');
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
fprintf(fid, '\n');

fprintf(fid, 'Files Generated:\n');
fprintf(fid, '- monthly_to_daily_results.mat: Main results structure\n');
fprintf(fid, '- complete_workspace.mat: All workspace variables\n');
fprintf(fid, '- daily_estimates.csv: Estimated daily values\n');
fprintf(fid, '- daily_estimates.png: Plot of daily estimates\n');
fprintf(fid, '- monthly_verification.png: Verification of monthly constraints\n');
fprintf(fid, '- error_metrics_comparison.png: Error metrics comparison\n');

fclose(fid);

fprintf('Process complete. All results saved to .mat files.\n');
fprintf('To load the results later, use: load(''monthly_to_daily_results.mat'')\n');
fprintf('To load the entire workspace, use: load(''complete_workspace.mat'')\n');
fprintf('A summary report has been saved to: summary_report.txt\n');
%% 5. Verification: Check if daily estimates aggregate correctly to monthly reports
fprintf('Verifying daily to monthly aggregation...\n');

% We already calculated monthly sums above, but let's enhance the verification
monthly_original = reports(:,3);
monthly_reconstructed = monthly_sums_non_int;

% Calculate discrepancy metrics
absolute_difference = monthly_original - monthly_reconstructed;
percentage_difference = (absolute_difference ./ monthly_original) * 100;

% Display summary statistics
fprintf('Verification Summary:\n');
fprintf('- Maximum Absolute Difference: %.2f\n', max(abs(absolute_difference)));
fprintf('- Mean Absolute Difference: %.2f\n', mean(abs(absolute_difference)));
fprintf('- Maximum Percentage Difference: %.2f%%\n', max(abs(percentage_difference)));
fprintf('- Mean Percentage Difference: %.2f%%\n', mean(abs(percentage_difference)));

% Count number of months with significant differences (e.g., >1%)
significant_diff_count = sum(abs(percentage_difference) > 1);
fprintf('- Number of months with >1%% difference: %d out of %d\n', significant_diff_count, length(monthly_original));

% Save detailed verification to a CSV file
verification_table = [
    (1:size(reports,1))', 
    reports(:,1), 
    reports(:,2), 
    monthly_original, 
    monthly_reconstructed, 
    absolute_difference, 
    percentage_difference
];
csvwrite('monthly_verification.csv', verification_table);

% Create a more detailed verification plot
figure;
subplot(2,1,1);
bar(absolute_difference);
title('Absolute Difference: Original - Reconstructed Monthly Values');
xlabel('Month Index');
ylabel('Difference');
grid on;

subplot(2,1,2);
bar(percentage_difference);
title('Percentage Difference: (Original - Reconstructed)/Original × 100%');
xlabel('Month Index');
ylabel('Percentage (%)');
grid on;
saveas(gcf, 'monthly_verification_detail.png');

% Add a scatter plot to visualize correlation between original and reconstructed
figure;
scatter(monthly_original, monthly_reconstructed, 50, 'filled');
hold on;
min_val = min(min(monthly_original), min(monthly_reconstructed));
max_val = max(max(monthly_original), max(monthly_reconstructed));
plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 1.5);
title('Original vs. Reconstructed Monthly Values');
xlabel('Original Monthly Reports');
ylabel('Sum of Daily Estimates');
grid on;
axis equal;
% Calculate correlation coefficient
corr_coef = corrcoef(monthly_original, monthly_reconstructed);
text(min_val + 0.1*(max_val-min_val), max_val - 0.1*(max_val-min_val), ...
     sprintf('Correlation: %.4f', corr_coef(1,2)), 'FontSize', 12);
saveas(gcf, 'monthly_correlation.png');

% Add these results to the structure
all_results.extended_verification = struct(...
    'monthly_original', monthly_original, ...
    'monthly_reconstructed', monthly_reconstructed, ...
    'absolute_difference', absolute_difference, ...
    'percentage_difference', percentage_difference, ...
    'max_abs_diff', max(abs(absolute_difference)), ...
    'mean_abs_diff', mean(abs(absolute_difference)), ...
    'max_pct_diff', max(abs(percentage_difference)), ...
    'mean_pct_diff', mean(abs(percentage_difference)), ...
    'significant_diff_count', significant_diff_count, ...
    'correlation', corr_coef(1,2) ...
);

% Update the summary report with the extended verification results
fid = fopen('summary_report.txt', 'a');
fprintf(fid, '\nExtended Verification Results:\n');
fprintf(fid, '----------------------------\n');
fprintf(fid, 'Maximum Absolute Difference: %.2f\n', max(abs(absolute_difference)));
fprintf(fid, 'Mean Absolute Difference: %.2f\n', mean(abs(absolute_difference)));
fprintf(fid, 'Maximum Percentage Difference: %.2f%%\n', max(abs(percentage_difference)));
fprintf(fid, 'Mean Percentage Difference: %.2f%%\n', mean(abs(percentage_difference)));
fprintf(fid, 'Number of months with >1%% difference: %d out of %d\n', significant_diff_count, length(monthly_original));
fprintf(fid, 'Correlation between original and reconstructed: %.4f\n', corr_coef(1,2));
fprintf(fid, '\nAdditional Files Generated:\n');
fprintf(fid, '- monthly_verification.csv: Detailed month-by-month verification\n');
fprintf(fid, '- monthly_verification_detail.png: Plots of differences\n');
fprintf(fid, '- monthly_correlation.png: Correlation plot\n');
fclose(fid);

% Re-save the updated all_results structure
save('monthly_to_daily_results.mat', 'all_results');

fprintf('Extended verification complete.\n');
fprintf('Additional verification results saved to monthly_verification.csv\n');

%% 6. Verification: Compare monthly reports with aggregated daily estimates
fprintf('Verifying that daily estimates aggregate back to monthly reports...\n');

% Create a table with monthly reports, sums of daily estimates, and differences
verif_table = zeros(size(reports,1), 4);

% For each month, calculate:
% 1. Original monthly report
% 2. Sum of daily estimates
% 3. Absolute difference
% 4. Percentage difference
for j = 1:size(reports,1)
    start_day = reports(j,1);
    end_day = reports(j,2);
    original_value = reports(j,3);
    
    % Sum the daily estimates for this month
    summed_value = sum(non_integer_daily(start_day:end_day));
    
    % Calculate differences
    abs_diff = original_value - summed_value;
    pct_diff = (abs_diff / original_value) * 100;
    
    % Store in table
    verif_table(j,1) = original_value;
    verif_table(j,2) = summed_value;
    verif_table(j,3) = abs_diff;
    verif_table(j,4) = pct_diff;
end

% Display summary statistics
total_original = sum(verif_table(:,1));
total_summed = sum(verif_table(:,2));
fprintf('Total of all monthly reports: %.2f\n', total_original);
fprintf('Total of all daily estimates: %.2f\n', total_summed);
fprintf('Overall difference: %.2f (%.4f%%)\n', total_original - total_summed,(total_original - total_summed)/total_original * 100);

% Calculate average absolute difference
mean_abs_diff = mean(abs(verif_table(:,3)));
mean_pct_diff = mean(abs(verif_table(:,4)));
fprintf('Mean absolute difference: %.2f (%.4f%%)\n', mean_abs_diff, mean_pct_diff);

% Count months with perfect alignment (difference < 0.01)
perfect_months = sum(abs(verif_table(:,3)) < 0.01);
fprintf('Months with perfect alignment (diff < 0.01): %d out of %d\n', perfect_months, size(reports,1));

% Save the verification table to CSV
headers = {'Month', 'Start_Day', 'End_Day', 'Original_Report', 'Sum_Daily_Estimates', 'Absolute_Difference', 'Percentage_Difference'};
month_indices = (1:size(reports,1))';
csv_table = [month_indices, reports(:,1:2), verif_table];
dlmwrite('monthly_verification_table.csv', csv_table, 'delimiter', ',');

% Also create a more readable text table
fid = fopen('monthly_comparison.txt', 'w');
fprintf(fid, 'Month\tStart\tEnd\tOriginal\tSum_Daily\tDifference\tPct_Diff\n');
fprintf(fid, '--------------------------------------------------------------------\n');
for j = 1:size(reports,1)
    fprintf(fid, '%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.4f%%\n',  j, reports(j,1), reports(j,2), verif_table(j,1), verif_table(j,2), verif_table(j,3), verif_table(j,4));
end
fprintf(fid, '--------------------------------------------------------------------\n');
fprintf(fid, 'Total:\t\t\t%.2f\t%.2f\t%.2f\t%.4f%%\n', total_original, total_summed, total_original - total_summed,  (total_original - total_summed)/total_original * 100);
fclose(fid);

% Display the table in a figure
figure;
t = uitable('Data', [month_indices, reports(:,1:2), verif_table], ...
           'ColumnName', headers, ...
           'RowName', {}, ...
           'Position', [20 20 750 400]);
title('Monthly Reports vs. Aggregated Daily Estimates');
saveas(gcf, 'monthly_comparison_table.png');

fprintf('Verification complete. Results saved to monthly_verification_table.csv and monthly_comparison.txt\n');
