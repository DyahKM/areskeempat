function [xhat_an, rmse, stop_an, L_stop, stop_rmse] = annihilating(A, y, x, events, ratio)
    % Get optimal filter length
    [rmse, stop_an, L_stop, stop_rmse] = stopLength(A, y, x, events);
    
    % Use the filter length corresponding to the given ratio
    L = L_stop(ratio);
    
    % Create annihilating filter matrix
    Xmat = zeros(length(x)-L+1, L);
    for s = 1:length(x)-L+1
        Xmat(s,:) = x(s:s+L-1);
    end
    
    Ymat = Xmat.'*Xmat;
    [U,Sigma,V] = svd(Ymat);
    h = U(:,end);
    c1 = [h(1); zeros(length(x)-L,1)];
    r1 = zeros(1,length(x));
    r1(1:L) = h.';
    Han = toeplitz(c1,r1);
    
    % Solve constrained optimization
    Aan = [A;Han];
    yan = [y; zeros(length(x)-L+1,1)];
    xhat_an = (pinv(Aan)*yan).';
    
    % Set RMSE to 0 since we don't need it
    rmse = 0;
end
    
    
