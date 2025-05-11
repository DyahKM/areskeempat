function variations = calculate_daily_variations(time_series)
% CALCULATE_DAILY_VARIATIONS Calculate daily variations in a time series
% Input:
%   time_series - Vector of daily values
% Output:
%   variations - Vector of daily variations (absolute differences)

% Calculate absolute differences between consecutive days
variations = abs(diff(time_series));

end 