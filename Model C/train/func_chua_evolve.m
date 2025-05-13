%clc; clear;
% New Parameters from the image
function [X,A,N,c3s]=func_chua_evolve(minx,maxx,miny,maxy,minz,maxz)
%clc; clear;

c1 = 15.6;
c2 = 1;
m0 = -8/7;
m1 = -5/7;

kappa = 0.3;    % Given parameter kappa
rng(6)
% Number of circuits
N = 6;
k=2;
a=31;
b=32;
%c3s=[50,35,33.8,33,25.58,40];
c3s = a + (b-a) * rand(N,1);
%c3s=[25,26,27,28,29];
% Varying parameters for each circuit, assuming they may differ as per original setup
% Defining them as arrays to accommodate different values for each circuit
c1s = c1 * ones(1, N); % If you want different c1 for each circuit, define them here
c2s = c2 * ones(1, N);
%c3s = c3 * ones(1, N);
m0s = m0 * ones(1, N);
m1s = m1 * ones(1, N);

% Initial conditions for all circuits (random small values near zero)

%x0 = 0.01 * randn(3*N, 1);
x0=0.01*ones(3*N, 1);

% Time span
tspan =0:0.001:500;

% Connection matrix
A = zeros(N, N);

    % 为每个节点添加k/2个连接
    for i = 1:N
        for j = 1:k/2
            % 计算连接的节点
            connect_to = mod(i+j-1, N) + 1;
            A(i, connect_to) = 1;
            A(connect_to, i) = 1;
        end
    end

% Solve the equations
[t, X] = ode45(@(t,x) coupled_chua_multi_adj(t, x, c1s, c2s, c3s, m0s, m1s, kappa, N, A), tspan, x0);

for m = 1:N
    X(:,3*m-2)=2*((X(:,3*m-2)-minx)/(maxx-minx))-1;
    X(:,3*m-1)=2*((X(:,3*m-1)-miny)/(maxy-miny))-1;
    X(:,3*m)=2*((X(:,3*m)-minz)/(maxz-minz))-1;
end

function dx = coupled_chua_multi_adj(t, x, c1s, c2s, c3s, m0s, m1s, kappa, N, A)
    dx = zeros(3*N, 1);
    % Calculate for each circuit
    for i = 1:N
        xi = x(3*i-2:3*i);
        g_xi = m1s(i)*xi(1) + (m0s(i) - m1s(i))*(abs(xi(1)+1) - abs(xi(1)-1))/2;
        
        % Compute coupling with other circuits
        coupling = 0;
        for j = 1:N
            if A(i, j) > 0
                xj = x(3*j-2);
                coupling = coupling + kappa * (xj - xi(1));
            end
        end
        
        dx1dt = c1s(i) * (xi(2) - xi(1) - g_xi) + coupling;
        dy1dt = c2s(i) * (xi(1) - xi(2) + xi(3));
        dz1dt = -c3s(i) * xi(2);
        
        dx(3*i-2:3*i) = [dx1dt; dy1dt; dz1dt];
    end
end
end
