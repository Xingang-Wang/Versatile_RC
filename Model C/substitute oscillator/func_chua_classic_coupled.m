%clc; clear;
% New Parameters from the image
function [X,A,N,c3s,minx,maxx,miny,maxy,minz,maxz]=func_chua_classic_coupled()
%clc; clear;

    N = 6;
    kappa = 0.5; % 耦合系数
    k = 2; % 每个节点的连接数量
    
    % 初始化参数并随机化（在合理范围内波动）
    rng(6);
    sigma_base = 10;
    rho_base = 35;
    beta_base = 8/3;
    
    sigmas = sigma_base*ones(1,N);% + 2 * rand(1, N) - 1; % [9, 11] 区间
    rhos = rho_base + 1 * rand(1, N) - 2; % [26, 30] 区间
    betas = beta_base*ones(1,N);%+ 0.5 * rand(1, N) - 0.25; % [7.75/3, 8.25/3] 区间
    c3s=rhos;
    % 初始条件（小的随机值）
    x0 = 0.01 * randn(3*N, 1);

    % 时间跨度
    tspan = 0:0.001:500;

    % 连接矩阵
    % A = zeros(N, N);
    % for i = 1:N
    %     for j = 1:k/2
    %         connect_to = mod(i+j-1, N) + 1;
    %         A(i, connect_to) = 1;
    %         A(connect_to, i) = 1;
    %     end
    % end
    % 创建一个N x N的零矩阵
% A = zeros(N, N);
% %设置第1个节点为中心节点，与其他所有节点相连
% A(1, 2:N) = 1;
% A(2:N, 1) = 1;
% Number of nodes
% Number of nodes
% Number of nodes
% Number of nodes
N = 6;

% Degrees of each node
d = [2, 3, 4, 2, 3, 4];

% Initialize the adjacency matrix
A = zeros(N);

% Create an array to keep track of remaining degree for each node
remaining_degrees = d;

% Generate the network
for i = 1:N
    % Connect current node to nodes with remaining degrees
    for j = i+1:N
        if remaining_degrees(i) > 0 && remaining_degrees(j) > 0
            A(i, j) = 1;
            A(j, i) = 1;
            remaining_degrees(i) = remaining_degrees(i) - 1;
            remaining_degrees(j) = remaining_degrees(j) - 1;
        end
    end
end

% Display the adjacency matrix
disp('Adjacency matrix:');
disp(A);
edge_function=[];

for i=1:size(A,1)
    for j=1:size(A,2)
        if A(i,j)==1
        edge_function=[edge_function;i,j];
        end
    end
end
% 仅提供文件名，保存到当前工作目录
filename1 = 'edge_loz_d_N_6.csv';
writematrix(edge_function, filename1);
D=sum(A,1);

    % 解决微分方程
    [t, X] = ode45(@(t,x) coupled_lorenz_multi_adj(t, x, sigmas, betas, rhos, kappa, N, A), tspan, x0);
    
    % 归一化数据
    minx = min(X(:,1));
    maxx = max(X(:,1));
    miny = min(X(:,2));
    maxy = max(X(:,2));
    minz = min(X(:,3));
    maxz = max(X(:,3));
    for m = 1:N
        X(:,3*m-2) = 2*((X(:,3*m-2)-minx)/(maxx-minx))-1;
        X(:,3*m-1) = 2*((X(:,3*m-1)-miny)/(maxy-miny))-1;
        X(:,3*m) = 2*((X(:,3*m)-minz)/(maxz-minz))-1;
    end
    
    % 绘图
    for i=1:6
    figure
    plot(X(:,3*i-2),'b.-');
    hold on;
    end
    function dx = coupled_lorenz_multi_adj(~, x, sigmas, betas, rhos, kappa, N, A)
        dx = zeros(3*N, 1);
        for i = 1:N
            xi = x(3*i-2:3*i);
            
            % 洛伦兹系统方程
            dx1dt = sigmas(i) * (xi(2) - xi(1));
            dy1dt = xi(1) * (rhos(i) - xi(3)) - xi(2);
            dz1dt = xi(1) * xi(2) - betas(i) * xi(3);
            
            % 耦合项
            coupling = 0;
            for j = 1:N
                if A(i, j) > 0
                    xj = x(3*j-2);
                    coupling = coupling + kappa * (xj - xi(1));
                end
            end
            dx(3*i-2:3*i) = [dx1dt + coupling; dy1dt; dz1dt];
        end
    end
end
