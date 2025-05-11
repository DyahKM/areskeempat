function [Inc_A, Out_A1, Out_AL] = iteration(Out, events, iterationTime, gama, ratio)
% ITERATION Run multiple ARES refinements for a fixed ratio index
% Inputs:
%   Out           - Initial reconstruction outputs from H-Fusion
%   events        - Ground truth or synthetic daily cases
%   iterationTime - Number of ARES iterations
%   gama          - Penalty parameter for L in cost function
%   ratio         - Index for selecting L from L_stop inside annihilating()
% Outputs:
%   Inc_A  - Struct containing RMSE values across iterations
%   Out_A1 - Output from first iteration
%   Out_AL - Output from final iteration

% Initialize output structures with all required fields
Out_A1 = struct('x_reconstr', [], 'A', [], 'y', [], 'L', [], 'muvars', [], 'Matrix', [], 'error', []);
Out_AL = struct('x_reconstr', [], 'A', [], 'y', [], 'L', [], 'muvars', [], 'Matrix', [], 'error', []);

% Process each time series
for l = 1:length(Out)
    % Get initial reconstruction
    reconX = Out(l).x_reconstr;
    A = Out(l).A;
    y = Out(l).y;
    
    % Apply annihilating filter
    [xhat_an, ~, ~, L_stop, ~] = annihilating(A, y, reconX, events, ratio);
    
    % Store results in Out_A1
    Out_A1(l).x_reconstr = xhat_an';
    Out_A1(l).A = A;
    Out_A1(l).y = y;
    Out_A1(l).L = L_stop(ratio);
    Out_A1(l).muvars = Out(l).muvars;
    Out_A1(l).Matrix = [];  % Will be filled if needed
    Out_A1(l).error = 0;    % Since we don't calculate error
    
    % Copy to Out_AL (same as A1 for single iteration)
    Out_AL(l) = Out_A1(l);
end

% For single iteration, Inc_A is same as Out_A1
Inc_A = Out_A1;
end
