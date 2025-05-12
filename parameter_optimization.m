function [best_params, best_score] = parameter_optimization(reports, num_days)
    % PARAMETER_OPTIMIZATION Tests different parameter combinations to find optimal values
    % Inputs:
    %   reports - Matrix containing monthly report data
    %   num_days - Total number of days in the dataset
    % Outputs:
    %   best_params - Structure containing the best parameter combination
    %   best_score - Score of the best parameter combination
    
    % Define parameter ranges to test
    lambda_range = [0.1, 0.5, 1, 2, 5, 10, 20];
    alpha_range = [0.1, 0.3, 0.5, 0.7, 0.9];
    gama_range = [0.1, 0.3, 0.5, 0.7, 0.9];
    iteration_range = [2, 4, 6, 8];
    stop_ratio_range = [4, 6, 8, 10, 12];
    
    % Initialize variables to store best results
    best_score = inf;
    best_params = struct();
    results = struct();
    
    % Create progress tracking
    total_combinations = length(lambda_range) * length(alpha_range) * length(gama_range) * ...
                        length(iteration_range) * length(stop_ratio_range);
    
    % Create figure for real-time visualization
    figure('Name', 'Parameter Optimization Progress', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    
    % Subplot for progress
    subplot(2,2,1);
    progress_bar = bar(0);
    title('Overall Progress');
    ylabel('Percentage Complete');
    ylim([0 100]);
    
    % Subplot for best score history
    subplot(2,2,2);
    best_scores = [];
    score_plot = plot(best_scores, 'b-', 'LineWidth', 2);
    title('Best Score History');
    xlabel('Combination Number');
    ylabel('Score');
    grid on;
    
    % Subplot for current parameters
    subplot(2,2,3);
    param_text = text(0.1, 0.5, 'Initializing...', 'FontSize', 10);
    title('Current Parameters');
    axis off;
    
    % Subplot for best parameters
    subplot(2,2,4);
    best_param_text = text(0.1, 0.5, 'No best parameters yet', 'FontSize', 10);
    title('Best Parameters Found');
    axis off;
    
    % Create log file
    log_fid = fopen('optimization_log.txt', 'w');
    fprintf(log_fid, 'Parameter Optimization Log\n');
    fprintf(log_fid, '=======================\n\n');
    fprintf(log_fid, 'Total combinations to test: %d\n\n', total_combinations);
    
    current_combination = 0;
    
    % Test all combinations
    for lambda = lambda_range
        for alpha = alpha_range
            for gama = gama_range
                for iterationTime = iteration_range
                    for stop_ratio = stop_ratio_range
                        current_combination = current_combination + 1;
                        
                        % Update progress visualization
                        progress = 100 * current_combination / total_combinations;
                        set(progress_bar, 'YData', progress);
                        drawnow;
                        
                        % Update current parameters display
                        current_params = sprintf(['Current Parameters:\n' ...
                            'Lambda: %.1f\n' ...
                            'Alpha: %.1f\n' ...
                            'Gama: %.1f\n' ...
                            'Iteration Time: %d\n' ...
                            'Stop Ratio: %d'], ...
                            lambda, alpha, gama, iterationTime, stop_ratio);
                        set(param_text, 'String', current_params);
                        
                        % Log current combination
                        fprintf(log_fid, '\nTesting combination %d/%d (%.1f%%)\n', ...
                            current_combination, total_combinations, progress);
                        fprintf(log_fid, 'Parameters: lambda=%.1f, alpha=%.1f, gama=%.1f, iter=%d, stop=%d\n', ...
                            lambda, alpha, gama, iterationTime, stop_ratio);
                        
                        try
                            % Get constraint matrix
                            [A, y] = rep_constraint_equations_full(reports, num_days);
                            
                            % Create initial reconstruction
                            temp_events = zeros(num_days, 1);
                            [recon_events, recon_error, ~, M] = sp_reconstruct(A, y, lambda, temp_events, alpha);
                            
                            % Create Out structure
                            Out = struct();
                            Out(1).muvars = [0, 0];
                            Out(1).A = A;
                            Out(1).y = y;
                            Out(1).x_reconstr = recon_events(:, 1, 1);
                            Out(1).x_error = recon_error;
                            Out(1).Matrix = M;
                            
                            % Run ARES filter
                            [~, ~, Out_AL] = iteration(Out, num_days, iterationTime, gama, stop_ratio);
                            
                            % Get daily estimates
                            daily_estimates = Out_AL(1).x_reconstr;
                            
                            % Apply integer constraints
                            integer_daily = enforce_integer_constraints(daily_estimates, reports, true);
                            
                            % Calculate error metrics
                            error_metrics = calculate_error_metrics(integer_daily, reports);
                            
                            % Calculate monthly sums for verification
                            monthly_sums = zeros(size(reports,1), 1);
                            for j = 1:size(reports,1)
                                monthly_sums(j) = sum(integer_daily(reports(j,1):reports(j,2)));
                            end
                            
                            % Calculate constraint validation
                            constraint_error = reports(:,3) - monthly_sums;
                            relative_error = abs(constraint_error) ./ reports(:,3) * 100;
                            
                            % Calculate overall score
                            score = mean(relative_error) + ...
                                   error_metrics.monthly_avg + ...
                                   error_metrics.neighboring + ...
                                   error_metrics.weekly + ...
                                   error_metrics.monthly;
                            
                            % Log metrics
                            fprintf(log_fid, 'Metrics:\n');
                            fprintf(log_fid, '- Monthly Average Error: %.4f\n', error_metrics.monthly_avg);
                            fprintf(log_fid, '- Neighboring Days Error: %.4f\n', error_metrics.neighboring);
                            fprintf(log_fid, '- Weekly Pattern Error: %.4f\n', error_metrics.weekly);
                            fprintf(log_fid, '- Monthly Pattern Error: %.4f\n', error_metrics.monthly);
                            fprintf(log_fid, '- Mean Relative Error: %.4f%%\n', mean(relative_error));
                            fprintf(log_fid, '- Overall Score: %.4f\n', score);
                            
                            % Store results
                            param_key = sprintf('lambda_%.1f_alpha_%.1f_gama_%.1f_iter_%d_stop_%d', ...
                                lambda, alpha, gama, iterationTime, stop_ratio);
                            results.(param_key) = struct(...
                                'params', struct('lambda', lambda, 'alpha', alpha, ...
                                               'gama', gama, 'iterationTime', iterationTime, ...
                                               'stop_ratio', stop_ratio), ...
                                'score', score, ...
                                'error_metrics', error_metrics, ...
                                'constraint_error', mean(relative_error));
                            
                            % Update best parameters if better score found
                            if score < best_score
                                best_score = score;
                                best_params = struct(...
                                    'lambda', lambda, ...
                                    'alpha', alpha, ...
                                    'gama', gama, ...
                                    'iterationTime', iterationTime, ...
                                    'stop_ratio', stop_ratio);
                                
                                % Update best parameters display
                                best_params_str = sprintf(['Best Parameters:\n' ...
                                    'Lambda: %.1f\n' ...
                                    'Alpha: %.1f\n' ...
                                    'Gama: %.1f\n' ...
                                    'Iteration Time: %d\n' ...
                                    'Stop Ratio: %d\n' ...
                                    'Score: %.4f'], ...
                                    lambda, alpha, gama, iterationTime, stop_ratio, score);
                                set(best_param_text, 'String', best_params_str);
                                
                                % Update score history plot
                                best_scores = [best_scores; score];
                                set(score_plot, 'XData', 1:length(best_scores), 'YData', best_scores);
                                ylim([min(best_scores)*0.9, max(best_scores)*1.1]);
                                
                                % Log new best parameters
                                fprintf(log_fid, '\n*** NEW BEST PARAMETERS FOUND ***\n');
                                fprintf(log_fid, 'Score: %.4f\n', score);
                                fprintf(log_fid, 'Parameters: lambda=%.1f, alpha=%.1f, gama=%.1f, iter=%d, stop=%d\n', ...
                                    lambda, alpha, gama, iterationTime, stop_ratio);
                            end
                            
                        catch e
                            fprintf(log_fid, 'Error: %s\n', e.message);
                            continue;
                        end
                    end
                end
            end
        end
    end
    
    % Save all results
    save('parameter_optimization_results.mat', 'results', 'best_params', 'best_score');
    
    % Create summary report
    fid = fopen('parameter_optimization_summary.txt', 'w');
    fprintf(fid, 'Parameter Optimization Summary\n');
    fprintf(fid, '===========================\n\n');
    fprintf(fid, 'Best Parameters Found:\n');
    fprintf(fid, 'Lambda: %.1f\n', best_params.lambda);
    fprintf(fid, 'Alpha: %.1f\n', best_params.alpha);
    fprintf(fid, 'Gama: %.1f\n', best_params.gama);
    fprintf(fid, 'Iteration Time: %d\n', best_params.iterationTime);
    fprintf(fid, 'Stop Ratio: %d\n', best_params.stop_ratio);
    fprintf(fid, '\nBest Score: %.4f\n', best_score);
    
    % Sort and display top 10 parameter combinations
    fprintf(fid, '\nTop 10 Parameter Combinations:\n');
    fprintf(fid, '----------------------------\n');
    
    % Convert results to array for sorting
    param_names = fieldnames(results);
    scores = zeros(length(param_names), 1);
    for i = 1:length(param_names)
        scores(i) = results.(param_names{i}).score;
    end
    
    [sorted_scores, idx] = sort(scores);
    for i = 1:min(10, length(sorted_scores))
        param = results.(param_names{idx(i)}).params;
        fprintf(fid, '%d. Score: %.4f\n', i, sorted_scores(i));
        fprintf(fid, '   Lambda: %.1f, Alpha: %.1f, Gama: %.1f\n', ...
            param.lambda, param.alpha, param.gama);
        fprintf(fid, '   Iteration Time: %d, Stop Ratio: %d\n\n', ...
            param.iterationTime, param.stop_ratio);
    end
    
    fclose(fid);
    fclose(log_fid);
    
    fprintf('\nOptimization complete!\n');
    fprintf('Results saved to parameter_optimization_results.mat\n');
    fprintf('Summary saved to parameter_optimization_summary.txt\n');
    fprintf('Detailed log saved to optimization_log.txt\n');
    
    % Keep the figure open
    uiwait(gcf);
end 