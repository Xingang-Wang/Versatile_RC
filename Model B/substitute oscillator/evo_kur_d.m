%% 2024.3.7
%clc;clear;
%Parameters

function [theta_sin,theta_cos,omega,A,d]=evo_kur_d()
N = 50; % Number of nodes
rng(2);
d= randi([2,4],50,1); % Desired degree
gamma = 5; % Assortativity parameter
delta=0.8;
[A,omega]=func_connectivity_many(N,gamma,delta,d);
D=mean(d);
edge_function=[];

for i=1:size(A,1)
    for j=1:size(A,2)
        if A(i,j)==1
        edge_function=[edge_function;i,j];
        end
    end
end
% 仅提供文件名，保存到当前工作目录
filename1 = 'edge_kur_d_N_50.csv';

% 写入 .csv 文件到当前工作目录
writematrix(edge_function, filename1);
%%

% Parameters
K = 0.4; % Coupling strength
T =2000; % Total simulation time
dt = 0.01; % Time step
timesteps = ceil(T / dt); % Number of time steps

% Initialize oscillator phases
theta = rand(N, 1) * 2 * pi - pi; % Random initial phases

% Initialize order parameter
R = zeros(timesteps, 1);
theta_sin = zeros(N,timesteps);
theta_cos = zeros(N,timesteps);
% Simulation loop using Runge-Kutta method
for t = 1:timesteps
    % Calculate k1
    k1 = dt * kuramoto_model(theta, omega, A, K);
    
    % Calculate k2
    k2 = dt * kuramoto_model(theta + k1 / 2, omega, A, K);
    
    % Calculate k3
    k3 = dt * kuramoto_model(theta + k2 / 2, omega, A, K);
    
    % Calculate k4
    k4 = dt * kuramoto_model(theta + k3, omega, A, K);
    
    % Update phases using RK4 formula
    theta = theta + (k1 + 2*k2 + 2*k3 + k4) / 6;
    
    % Wrap angles to [-pi, pi)
    theta = mod(theta + pi, 2*pi) - pi;
    theta_sin(:,t)=sin(theta);
    theta_cos(:,t)=cos(theta);
    % Calculate order parameter
    R(t) = abs(sum(A * exp(1i * theta))) / (N * D);
end
%save Kur_sin_cos_A_rng2.mat theta_sin theta_cos omega A d
% Plot the evolution of the order parameter
figure
t = linspace(0, T, length(R));
plot(t, R);
xlabel('Time');
ylabel('Order parameter R(t)');
title('Evolution of the order parameter');
figure
plot(t(1:5000),theta_sin(3,1:5000));


% Function to compute the Kuramoto model equation

