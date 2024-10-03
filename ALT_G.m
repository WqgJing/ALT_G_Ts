function [X, Q_X] = ALT_G(Q_0,  k, rho_s, C_s, theta_i, rho_i, L,dt)


    alpha = k / (rho_s * C_s);          % Thermal diffusivity
    theta_rho_L = theta_i * rho_i * L;  % Constant for heat flux calculation
    coeff = (2 * k) / (rho_s * C_s);    % Constant for Q_X

    % Time discretization
    T_total = numel(Q_0) * dt;  % Total simulation time
    N = numel(Q_0);                     
    t = linspace(0, T_total, N+1);     

    % Initialize variables
    X = zeros(1, N+1);                  % Position X(t)
    X(1) = 0.01;                        % Initial position
    Q_X = zeros(1, N+1);                % Heat flux at X(t)

    % Simulation loop
    for n = 2:N+1
        % Initial guess for X(n)
        X_n_old = X(n-1);
        iter = 0;
        max_iter = 20;
        tol = 1e-6;
        converged = false;
        
        while ~converged && iter < max_iter
            I_n = 0;
            for j = 1:n-1
                t_diff = t(n) - t(j);
                if t_diff == 0
                    continue; 
                end
                
                S_nj_X = (1 / sqrt(pi * 4*alpha * t_diff)) * ...
                   ((X_n_old - X(j))./(2 * alpha * t_diff) * (exp(-((X_n_old - X(j))^2) / (4 * alpha * t_diff)) + ...
                    (X_n_old + X(j))./(2 * alpha * t_diff) * exp(-((X_n_old + X(j))^2) / (4 * alpha * t_diff))));
                 
                S_nj_0 = (2 / sqrt(pi * 4*alpha * t_diff)) * ...
                    (X_n_old )./(2 * alpha * t_diff) * exp(-((X_n_old)^2) / (4 * alpha * t_diff));
                
                if j == 1
                    Q_X(j) = theta_rho_L * (X(j) - X(j)) / dt; % First term zero
                else
                    Q_X(j) = theta_rho_L * (X(j) - X(j-1)) / dt;
                end
                
                I_n = I_n + (Q_X(j) * S_nj_X - Q_0(j) * S_nj_0) * dt;
            end
            
            Q_X_n = -coeff * I_n;
            
            X_n_new = X(n-1) + (dt / theta_rho_L) * Q_X_n;
            
            if abs(X_n_new - X_n_old) < tol
                converged = true;
            else
                X_n_old = X_n_new;
                iter = iter + 1;
            end
        end
        
        X(n) = X_n_new;
        Q_X(n) = Q_X_n;
    end
    X=X-0.01;
    X(1)=[];
    Q_X(1)=[];
end
