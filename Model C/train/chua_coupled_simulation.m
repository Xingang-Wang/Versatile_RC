%function chua_coupled_simulation()
clc;clear;
    % 电路参数
    C3=10.8e-06;
    C4=88.2e-06;
    C5=10.9e-06;
    C6=94.6e-06;
    R2=1.6e3;
    R3=1.6e3;
    R4=10;
    R5=10;
    %R6=13.72e3;
    R6=0;
    L1=26.7e-3;
    L2=26.5e-3;

     % C3=10*10^-9;
     % C5=10*10^-9;
     % C6=100*10^-9;
     % C4=100*10^-9;
     % L1=26*10^-3;
     % L2=26*10^-3;
     % R2=1600;
     % R3=1600;
     % R4=10;R5=10;
     % %R6=150000;


     tspan = [0:1*10^-8:1*10^-3]; % 无量纲时间范围
    % 初始条件：两个电路的[v'_1, v'_2, i'_1L]
    init_cond = rand(6,1);
    %init_cond = zeros(6,1);
    %init_cond = [-0.5;0;0.2;-0.5;0;0.2];
    % ODE求解器
   % 设置误差容许范围
options = odeset('RelTol',1e-5,'AbsTol',[1e-6 1e-6 1e-6 1e-6 1e-6 1e-6]);

% 使用ode45求解
[T, Y] = ode45(@(t, y) coupled_chua_ode(t, y, C3, C4, C5, C6, R2, R3, R4, R5, R6, L1, L2), tspan, init_cond, options);
    
    % 绘图
    figure;
    %subplot(2,1,1);
    plot( Y(:,1), 'r-');
    hold on
    plot(Y(:,4),'b-');
    title('v''_1 (C1 Voltage) for Both Circuits');
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    legend('Circuit 1', 'Circuit 2');
    
    


function dydt = coupled_chua_ode(t, y,C3,C4,C5,C6,R2,R3,R4,R5,R6,L1,L2)
    % y = [v'_1_1, v'_2_1, i'_L_1, v'_1_2, v'_2_2, i'_L_2]
    y=rand(6,1);
    % 微分方程组，包括耦合
    dv1dt_1 = 1/C3*(1/R2*(y(1)-y(2))-chua_diode(y(1))+1/R6*(y(4)-y(1)));
    dv2dt_1 = 1/C4*(1/R2*(y(2)-y(1))+y(3));
    diLdt_1 = 1/L1*(-y(2)-R4*y(3));
    
    dv1dt_2 = 1/C5*(1/R3*(y(4)-y(5))-chua_diode(y(4))+1/R6*(y(1)-y(4)));
    dv2dt_2 = 1/C6*(1/R3*(y(5)-y(4))+y(6));
    diLdt_2 =  1/L2*(-y(5)-R5*y(6));
    
    dydt = [dv1dt_1; dv2dt_1; diLdt_1; dv1dt_2; dv2dt_2; diLdt_2];
end

% Chua二极管的非线性函数定义如上例
function g = chua_diode(v)
Bp=1.7;
m0=-0.41*10^-3;
m1=-0.76*10^-3;
g=m0*v+0.5*(m1-m0)*(abs(v+Bp)-abs(v-Bp));
%y=0.8041
% E=1;
% Ga=-0.6955e-3;
% Gb=-0.4347e-3;
%     % Chua二极管的非线性函数，根据实际情况调整参数
%v=0.8041;
%g=Gb*v+0.5*(Ga-Gb)*(abs(0.8041+E)-abs(0.8041-E));
end

