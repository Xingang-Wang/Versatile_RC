clc; clear;

c1 = 15.6;
c2 = 1;
m0 = -8/7;
m1 = -5/7;
kappa = 0.2;    % Given parameter kappa
rng(6)
N = 6;
k=2;
a=35;
b=36;
%c3s=[50,35,33.8,33,25.58,40];
c3s = a + (b-a) * rand(N,1);
c1s = c1 * ones(1, N); 
c2s = c2 * ones(1, N);
m0s = m0 * ones(1, N);
m1s = m1 * ones(1, N);

x0 = 0.01 * randn(3*N, 1);

tspan = [0 :0.02: 500];
h = 0.01; % Step size for RK4
% 
% A = zeros(N, N);
% for i = 1:N
%     for j = 1:k/2
%         connect_to = mod(i+j-1, N) + 1;
%         A(i, connect_to) = 1;
%         A(connect_to, i) = 1;
%     end
% end
% 创建一个N x N的零矩阵
A = zeros(N, N);

% 设置第1个节点为中心节点，与其他所有节点相连
A(1, 2:N) = 1;
A(2:N, 1) = 1;
% Preallocate space for solution and order parameters
num_steps = length(tspan);
X = zeros(num_steps, 3*N);
X(1, :) = x0';
order_params = zeros(num_steps, 1); % Store the order parameter

for idx = 1:num_steps-1
    t = tspan(idx);
    y = X(idx, :)';
    
    k1 = h * coupled_chua_multi_adj(t, y, c1s, c2s, c3s, m0s, m1s, kappa, N, A);
    k2 = h * coupled_chua_multi_adj(t + h/2, y + k1/2, c1s, c2s, c3s, m0s, m1s, kappa, N, A);
    k3 = h * coupled_chua_multi_adj(t + h/2, y + k2/2, c1s, c2s, c3s, m0s, m1s, kappa, N, A);
    k4 = h * coupled_chua_multi_adj(t + h, y + k3, c1s, c2s, c3s, m0s, m1s, kappa, N, A);
    
    y_next = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
    X(idx+1, :) = y_next';

    % Calculate average Euclidean distance for order parameter
    total_distance = 0;
    count = 0;
    for i = 1:N
        for j = i+1:N
            xi = y_next(3*i-2:3*i);
            xj = y_next(3*j-2:3*j);
            distance = norm(xi - xj);
            total_distance = total_distance + distance;
            count = count + 1;
        end
    end
    order_params(idx+1) = total_distance / count; % Store the average distance
end
  minx = min(X(:,1));
    maxx = max(X(:,1));
    miny = min(X(:,2));
    maxy = max(X(:,2));
    minz = min(X(:,3));
    maxz = max(X(:,3));
for m = 1:N
    X(:,3*m-2)=2*((X(:,3*m-2)-minx)/(maxx-minx))-1;
    X(:,3*m-1)=2*((X(:,3*m-1)-miny)/(maxy-miny))-1;
    X(:,3*m)=2*((X(:,3*m)-minz)/(maxz-minz))-1;
end
% Plotting the results
figure;
plot(tspan, order_params);
title('Average Euclidean Distance over Time');
xlabel('Time');
ylabel('Average Distance');

figure;
plot(tspan, X(:, 1:3:end));
title('State Variable x over Time for All Circuits');
xlabel('Time');
ylabel('State Variable x');

function dx = coupled_chua_multi_adj(t, x, c1s, c2s, c3s, m0s, m1s, kappa, N, A)
    dx = zeros(3*N, 1);
    for i = 1:N
        xi = x(3*i-2:3*i);
        g_xi = m1s(i)*xi(1) + (m0s(i) - m1s(i))*(abs(xi(1)+1) - abs(xi(1)-1))/2;
        
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
