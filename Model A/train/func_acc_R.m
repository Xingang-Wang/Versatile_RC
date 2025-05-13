%% 2024.3.7
%clc;clear;
function [acc_R,theta_sin]=func_acc_R(rng_num)
N = 50; % Number of nodes
d = 3; % Desired degree
gamma = 5; % Assortativity parameter
delta=0.8;
[A,omega]=func_connectivity_A(N,gamma,delta,d);
%%
rng(rng_num);
% Parameters
K =0.4; % Coupling strength
T =150; % Total simulation time
dt = 0.1; % Time step
timesteps = ceil(T / dt); % Number of time steps

% Initialize oscillator phases
 % Random initial phases
 theta = rand(N, 1) * 2 * pi - pi;
% Initialize order parameter
acc_R = zeros(timesteps, 1);
theta_sin = zeros(N,timesteps);
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
    % Calculate order parameter
    acc_R(t) = abs(sum(A * exp(1i * theta))) / (N * d);
end
% figure
% plot(acc_R,'r.-')