% Test different RMSE calculations using baselinemain8.m output
% This file analyzes the daily estimates using 4 different RMSE approaches

% Load the results from baselinemain8.m
load('monthly_to_daily_results.mat');
load('complete_workspace.mat');
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

% Get the daily estimates
daily_estimates = all_results.final_daily_estimates;

%% Create date range (from Jan 1, 2020)
start_date = datetime(2020,1,1);
dates = start_date + (0:length(daily_estimates)-1)';

% Create table with dates and estimates
daily_table = table(dates, daily_estimates, 'VariableNames', {'Date', 'Daily_Estimates'});

% Export to CSV
csv_filename = 'daily_estimates.csv';
writetable(daily_table, csv_filename);
fprintf('✓ Exported to CSV: %s\n', csv_filename);

% Export to Excel
excel_filename = 'daily_estimates.xlsx';
writetable(daily_table, excel_filename);
fprintf('✓ Exported to Excel: %s\n', excel_filename);

% Add monthly information to Excel
monthly_table = table();
for j = 1:size(reports,1)
    start_day = reports(j,1);
    end_day = reports(j,2);
    month_data = daily_estimates(start_day:end_day);
    monthly_table.Date(j) = dates(start_day);
    monthly_table.Monthly_Total(j) = reports(j,3);
    monthly_table.Sum_of_Daily(j) = sum(month_data);
    monthly_table.Difference(j) = reports(j,3) - sum(month_data);
    monthly_table.Percentage_Error(j) = abs(monthly_table.Difference(j)) / reports(j,3) * 100;
end

% Write monthly summary to second sheet of Excel
writetable(monthly_table, excel_filename, 'Sheet', 'Monthly Summary');
fprintf('✓ Added monthly summary to Excel file\n');

%% 1. RMSE based on Monthly Averages
% This measures how well daily estimates match the expected daily average for each month
fprintf('\n1. ERROR based on Monthly Averages:\n');
fprintf('This measures deviation from expected daily averages\n');

% Calculate average daily value for each month
monthly_avg = reports(:,3) ./ (reports(:,2) - reports(:,1) + 1);

% Calculate RMSE
monthly_avg_errors = zeros(size(reports,1), 1);
for j = 1:size(reports,1)
    % Get daily estimates for this month
    month_data = daily_estimates(reports(j,1):reports(j,2));
    % Calculate average of daily estimates
    daily_avg = mean(month_data);
    % Calculate squared error
    monthly_avg_errors(j) = (monthly_avg(j) - daily_avg)^2;
end
error_monthly_avg = sqrt(mean(monthly_avg_errors));
fprintf('error_monthly_avg: %.4f\n', error_monthly_avg);
fprintf('Interpretation: Lower value means daily estimates are closer to expected monthly averages\n');

%% 2. RMSE based on Neighboring Days
% This measures the smoothness of transitions between consecutive days
fprintf('\n2. ERROR based on Neighboring Days:\n');
fprintf('This measures smoothness between consecutive days\n');

% Calculate differences between consecutive days
day_differences = diff(daily_estimates);
% Calculate RMSE
error_neighboring = sqrt(mean(day_differences.^2));
fprintf('Error_neighboring: %.4f\n', error_neighboring);
fprintf('Interpretation: Lower value means smoother transitions between days\n');

%% 3. RMSE based on Weekly Patterns
% This measures the consistency of patterns within each week
fprintf('\n3. ERROR based on Weekly Patterns:\n');
fprintf('This measures consistency within weekly patterns\n');

% Calculate weekly variances
num_weeks = floor(length(daily_estimates)/7);
weekly_errors = zeros(num_weeks, 1);
for w = 1:num_weeks
    % Get data for this week
    week_data = daily_estimates((w-1)*7+1:w*7);
    % Calculate variance within week
    weekly_errors(w) = var(week_data);
end
error_weekly = sqrt(mean(weekly_errors));
fprintf('error_weekly: %.4f\n', error_weekly);
fprintf('Interpretation: Lower value means more consistent weekly patterns\n');

%% 4. RMSE based on Monthly Patterns
% This measures the consistency of patterns within each month
fprintf('\n4. ERROR based on Monthly Patterns:\n');
fprintf('This measures consistency within monthly patterns\n');

% Calculate monthly standard deviations
monthly_std = zeros(size(reports,1), 1);
for j = 1:size(reports,1)
    % Get data for this month
    month_data = daily_estimates(reports(j,1):reports(j,2));
    % Calculate standard deviation within month
    monthly_std(j) = std(month_data);
end
error_monthly = sqrt(mean(monthly_std));
fprintf('error_monthly: %.4f\n', error_monthly);
fprintf('Interpretation: Lower value means more consistent monthly patterns\n');

%% Save Results
% Create a structure to store all RMSE results
error_results = struct();
error_results.monthly_avg = error_monthly_avg;
error_results.neighboring = error_neighboring;
error_results.weekly = error_weekly;
error_results.monthly = error_monthly;

% Save results
save('error_analysis_results.mat', 'error_results');

% Create a summary table
fprintf('\nSummary of RMSE Analysis:\n');
fprintf('------------------------\n');
fprintf('1. Monthly Average RMSE: %.4f\n', error_monthly_avg);
fprintf('2. Neighboring Days RMSE: %.4f\n', error_neighboring);
fprintf('3. Weekly Pattern RMSE: %.4f\n', error_weekly);
fprintf('4. Monthly Pattern RMSE: %.4f\n', error_monthly);

% Plot results
figure;
rmse_values = [error_monthly_avg, error_neighboring, error_weekly, error_monthly];
bar(rmse_values);
title('ERROR Analysis Results');
xlabel('ERROR Type');
ylabel('ERROR Value');
set(gca, 'XTickLabel', {'Monthly Avg', 'Neighboring', 'Weekly', 'Monthly'});
grid on;
saveas(gcf, 'rmse_analysis_plot.png'); 