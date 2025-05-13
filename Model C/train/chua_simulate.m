clc;clear;
% Parameters
a = -1.0914;    % Given parameter a
b = 0.0579;    % Given parameter b
alpha1 = 5.4; % Given parameter alpha
alpha2 = 5.9;
beta = 6.0194;   % Given parameter beta
kappa = 0.3;     % Coupling strength
x0 = [0.01; 0.0; 0; 0.01; 0; 0]; % Initial conditions for both circuits
%x0=rand(6,1);1
% Time span
tspan = [0 200];

% Solve the equations
[t, X] = ode45(@(t,x) coupled_chua(t, x, a, b, alpha1,alpha2, beta, kappa), tspan, x0);

% Plot results
figure;
plot(t, X(:,1), t, X(:,4));
legend('x1', 'x2');
%title('Coupled Chua Circuits');
xlabel('Time');
ylabel('State Variables');



function dx = coupled_chua(t, x, a, b, alpha1,alpha2, beta, kappa)
    % x = [x1, y1, z1, x2, y2, z2]'

    % Chua's diode characteristic for both circuits
    h1 = b*x(1) + 0.5*(a - b)*(abs(x(1)+1) - abs(x(1)-1));
    h2 = b*x(4) + 0.5*(a - b)*(abs(x(4)+1) - abs(x(4)-1));
    
    % System equations
    dx1dt = alpha1 * (x(2) - x(1) - h1) + kappa * (x(4) - x(1));
    dy1dt = x(1) - x(2) + x(3);
    dz1dt = -beta * x(2);
    
    dx2dt = alpha2 * (x(5) - x(4) - h2) + kappa * (x(1) - x(4));
    dy2dt = x(4) - x(5) + x(6);
    dz2dt = -beta * x(5);
    
    dx = [dx1dt; dy1dt; dz1dt; dx2dt; dy2dt; dz2dt];
end


