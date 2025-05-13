clc; clear;
% Parameters
a = -1.0914;    % Given parameter a
b = 0.0579;     % Given parameter b
betas = 6.0194; % Given parameter beta
kappa = 0.1;    % Coupling strength

% Number of circuits
N = 5;

% Different alpha for each circuit
alphas = linspace(5.4, 5.9, N);

% Initial conditions for all circuits (random small values near zero)
x0 = 0.01 * randn(3*N, 1);

% Time span
tspan = [0 200];

% Connection matrix
A = zeros(N, N);
A(1, 2:N) = 1; % Central node connects to all other nodes
A(2:N, 1) = 1; % Symmetric connection

% Solve the equations
[t, X] = ode45(@(t,x) coupled_chua_multi_adj(t, x, a, b, alphas, betas, kappa, N, A), tspan, x0);

% Plot results
figure;
hold on;
for i = 1:N
    plot(t, X(:, 3*i-2)); % Plot x variable for each circuit
end
hold off;
legend(arrayfun(@(i) sprintf('x%d', i), 1:N, 'UniformOutput', false));
xlabel('Time');
ylabel('State Variables');

function dx = coupled_chua_multi_adj(t, x, a, b, alphas, betas, kappa, N, A)
    dx = zeros(3*N, 1);
    % Calculate for each circuit
    for i = 1:N
        xi = x(3*i-2:3*i);
        hi = b*xi(1) + 0.5*(a - b)*(abs(xi(1)+1) - abs(xi(1)-1));
        
        % Compute coupling with other circuits
        coupling = 0;
        for j = 1:N
            if A(i, j) > 0
                xj = x(3*j-2);
                coupling = coupling + kappa * (xj - xi(1));
            end
        end
        
        dx1dt = alphas(i) * (xi(2) - xi(1) - hi) + coupling;
        dy1dt = xi(1) - xi(2) + xi(3);
        dz1dt = -betas * xi(2);
        
        dx(3*i-2:3*i) = [dx1dt; dy1dt; dz1dt];
    end
end
