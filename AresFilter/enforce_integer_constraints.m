function constrained_daily = enforce_integer_constraints(daily_estimates, reports, enforce_integer)
% ENFORCE_INTEGER_CONSTRAINTS Apply constraints to daily estimates
% Inputs:
%   daily_estimates - Vector of daily estimates
%   reports        - Matrix of monthly reports [start_day, end_day, total]
%   enforce_integer- Boolean to enforce integer values
% Output:
%   constrained_daily - Daily estimates with constraints applied

% Initialize output
constrained_daily = daily_estimates;

% Process each month
for i = 1:size(reports, 1)
    start_day = reports(i, 1);
    end_day = reports(i, 2);
    monthly_total = reports(i, 3);
    
    % Get current month's estimates
    month_estimates = daily_estimates(start_day:end_day);
    
    if enforce_integer
        % For integer version:
        % 1. Round all values
        rounded = round(month_estimates);
        
        % 2. Calculate difference from monthly total
        current_sum = sum(rounded);
        diff = monthly_total - current_sum;
        
        % 3. Adjust values to match monthly total
        if diff ~= 0
            % Sort by decimal part to determine which values to adjust
            [~, idx] = sort(abs(month_estimates - rounded));
            
            % Add or subtract 1 from values with largest decimal parts
            for j = 1:abs(diff)
                if diff > 0
                    rounded(idx(j)) = rounded(idx(j)) + 1;
                else
                    rounded(idx(j)) = rounded(idx(j)) - 1;
                end
            end
        end
        
        % Store adjusted values
        constrained_daily(start_day:end_day) = rounded;
    else
        % For non-integer version:
        % Simply scale to match monthly total
        current_sum = sum(month_estimates);
        if current_sum ~= 0
            scale_factor = monthly_total / current_sum;
            constrained_daily(start_day:end_day) = month_estimates * scale_factor;
        end
    end
end

end 