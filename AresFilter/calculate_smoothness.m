function smoothness = calculate_smoothness(time_series)
% CALCULATE_SMOOTHNESS Calculate smoothness score of a time series
% Input:
%   time_series - Vector of daily values
% Output:
%   smoothness - Smoothness score (lower is smoother)

% Calculate first differences
diffs = diff(time_series);

% Calculate smoothness as the sum of squared differences
% Lower value means smoother series
smoothness = sum(diffs.^2) / length(diffs);

end 