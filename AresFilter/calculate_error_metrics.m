function error_metrics = calculate_error_metrics(daily_estimates, reports)
% CALCULATE_ERROR_METRICS Calculate various error metrics for daily estimates
% Inputs:
%   daily_estimates - Vector of daily estimates
%   reports        - Matrix of monthly reports [start_day, end_day, total]
% Output:
%   error_metrics - Structure containing various error metrics

% Initialize error metrics structure
error_metrics = struct('monthly_avg', 0, 'neighboring', 0, 'weekly', 0, 'monthly', 0);

% 1. Monthly Average Error
% Calculate average daily cases for each month
monthly_avgs = zeros(size(reports, 1), 1);
for i = 1:size(reports, 1)
    start_day = reports(i, 1);
    end_day = reports(i, 2);
    days_in_month = end_day - start_day + 1;
    monthly_avgs(i) = reports(i, 3) / days_in_month;
end

% Calculate error between daily estimates and monthly averages
monthly_avg_errors = zeros(size(reports, 1), 1);
for i = 1:size(reports, 1)
    start_day = reports(i, 1);
    end_day = reports(i, 2);
    month_estimates = daily_estimates(start_day:end_day);
    monthly_avg_errors(i) = mean(abs(month_estimates - monthly_avgs(i)));
end
error_metrics.monthly_avg = mean(monthly_avg_errors);

% 2. Neighboring Days Error
% Calculate average difference between consecutive days
daily_diffs = abs(diff(daily_estimates));
error_metrics.neighboring = mean(daily_diffs);

% 3. Weekly Pattern Error
% Calculate average difference between same weekdays
weekly_errors = zeros(7, 1);
for i = 1:7
    % Get all values for this weekday
    weekday_values = daily_estimates(i:7:end);
    % Calculate average difference from mean
    weekly_errors(i) = mean(abs(weekday_values - mean(weekday_values)));
end
error_metrics.weekly = mean(weekly_errors);

% 4. Monthly Pattern Error
% Calculate average difference between same days of different months
monthly_pattern_errors = zeros(31, 1);
for i = 1:31
    % Get all values for this day of month
    day_values = [];
    for j = 1:size(reports, 1)
        start_day = reports(j, 1);
        if i <= (reports(j, 2) - start_day + 1)
            day_values = [day_values; daily_estimates(start_day + i - 1)];
        end
    end
    if ~isempty(day_values)
        monthly_pattern_errors(i) = mean(abs(day_values - mean(day_values)));
    end
end
error_metrics.monthly = mean(monthly_pattern_errors(monthly_pattern_errors > 0));

end 