%% 2024.3.9
clc;clear;
load opt_data_kur_index_10_new.mat 
load min_rng_set_kur_index_10_new.mat
%d=3;
%%
 result = getfield ( opt_trials,'Fval');
 param= getfield ( opt_trials,'X');
 [sort_result,result_num]=sort(result);
 sort_param=param(result_num,:);
 opt_result=sort_param(1,:);
 sort_rng=min_rng_set(result_num);
 opt_rng=sort_rng(1);  
 rng(opt_rng);
 %%
 d=3;
max_d=3;
drive_num=50;
data_num=train_num;
testnum=train_num;
test_num=train_num;
data_len=ones(1,data_num)*train_length;
 %%
resSize =500; % size of the reservoir nodes;  
initLen = 100;
TrainLen=sum(data_len)-1;
test_Len =1000;
testLen=1000;
inSize = 2*(d+1)+1; 
outSize = 2;
nonliner_num=2;
X=[];
Yt = train_output_data(1:outSize,2:TrainLen+1);% run the reservoir with the data and collect X
%%

inSize=2*(d+1)+1;
W_in_a = opt_result(2);
Win= (2.0*rand(resSize,inSize)-1.0)*W_in_a;


%%
for i=1:train_num
    indata=train_data{i};
    [X1,W,reg]=func_get_X(inSize,opt_result,opt_rng,indata,Win,train_length);
    X=[X,X1];
end
X=X(:,1:end-1);
%% 
Data_len=[0,data_len];
for i=0:data_num-1
   trainLen=sum(Data_len(1:i+1));
   X(:,trainLen+1-i*initLen:trainLen+initLen-i*initLen)=[];
   Yt(:,trainLen+1-i*initLen:trainLen+initLen-i*initLen)=[];
end
rank=randperm( size(X,2) );  
X=X(:, rank); 
Yt=Yt(:, rank); 
X_T = X';
Wout = Yt*X_T / (X*X_T + reg*eye(nonliner_num*resSize+1));
a = opt_result(3);
 x1=2*rand(resSize,1)-1;
%%
%Parameters
N = 50; % Number of nodes
d = 3; % Desired degree
gamma = 5; % Assortativity parameter
delta=0.8;
A_rng=2;
[A,omega]=func_connectivity_A(N,gamma,delta,d);
omega_RC=omega;
% if min(omega) < 0
%     omega_RC = omega + abs(min(omega));
% else
%    omega_RC = omega; % 如果所有元素都是非负的，则不需要调整
% end
% %omega_RC=omega_RC*2;
% omega_RC=omega*3;% 
%%

% Parameters
K = 0.4; % Coupling strength
T =150; % Total simulation time
dt = 0.1; % Time step
timesteps = ceil(T / dt); % Number of time steps


%%
edge=[];
for m=1:size(A,1)
    for n=1:size(A,2)
        if A(m,n)==1
        edge=[edge;m,n,omega(m),omega(n)];
        end
    end
end

%%
rep_num=1;
rep_index =19;
rep_all_indexs=zeros(rep_num,d+1);
for i=1:rep_num
    rep_all_indexs(i,1)=rep_index(i);
    rep_all_indexs(i,2:end)=edge(edge(:,1)==rep_index(i),2);
end

%%
rng_num=9;
rng(rng_num);
test_data=zeros(2*(d+1)+1,1);
theta = rand(N, 1) * 2 * pi - pi; % Random initial phases;

% Simulation loop using Runge-Kutta method
for t = 1:timesteps
    theta = mod(theta + pi, 2*pi) - pi;
    theta_sin=sin(theta);
    theta_cos=cos(theta);
    %%
    for i=1:d+1
     test_data(2*i-1)=theta_sin(rep_all_indexs(i));
     test_data(2*i)=theta_cos(rep_all_indexs(i));
    end
     test_data(2*(d+1)+1)=omega_RC(rep_all_indexs(1,1));
     u=test_data;
     x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
     y = Wout*[1;x1;x1.^2;];
     RC_rep_theta=atan2(y(1), y(2));
     Y1(t)=y(1);
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
    
    %%
   
    
    if t>drive_num
        theta(rep_index)=RC_rep_theta;
    end
    % Calculate order parameter
    R(t) = abs(sum(A * exp(1i * theta))) / (N * d);
end
color_f1= addcolorplus(16);
color_f3= addcolorplus(3);
color_f2= addcolorplus(185);
%color_f2= addcolorplus(73);
color_f4= addcolorplus(215);
[acc_R,theta_sin]=func_acc_R(rng_num);
epsilon_percent=0.15;
tau=50;
%N_precise = func_calculate_precise_prediction_periods(acc_R, R, epsilon_percent, tau)

% Plot the evolution of the order parameter
t = linspace(0, T, timesteps)*0.115;
TightPlot.ColumeNumber = 2;     % 子图行数
TightPlot.RowNumber = 1;    % 子图列数
TightPlot.GapW = 0.2;  % 子图之间的左右间距
TightPlot.GapH = 0.0;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.25;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.12;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.18;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.015;  % 子图与图片右方的间距

%% PLOT
%figure(1);  % 声明Figure
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
figure
axes(p(1));
t_offset = t - t(drive_num+1);  % 偏移 t，使虚线位置为 0

% 重新绘制图像
% plot(t_offset(1:drive_num),theta_sin(rep_index,1:drive_num),'Marker','o', 'MarkerSize', 2, ...
%    'MarkerEdgeColor', color_f1,'MarkerFaceColor', color_f1);
% 
% hold on;
plot(t_offset(drive_num+1:end),theta_sin(rep_index,drive_num+1:timesteps),...
   'Linestyle','-','linewidth',3,'Color',color_f3);
hold on;
plot(t_offset(drive_num+1:end),Y1(drive_num+1:timesteps),...
   'Linestyle','-.','linewidth',1.5,'Color',color_f2);
set(gca,'FontName','Times New Roman','FontSize',13);

ylim([-1.05 1.05]);
%xlabel('$ \Lambda $','interpreter','latex');
%ylabel('$ sin(\theta)  $','interpreter','latex'
ylabel('\boldmath$ sin(\theta)  $','interpreter','latex','FontSize', 19);
% 调整 xlim，以虚线位置为0，范围重新设置为 [0, 原始t的最大值 - 虚线位置对应的t值]
%xlim([min(t_offset)-0.1 max(t_offset)+0.1-4.5]);
xlim([0 max(t_offset)+0.1-4.5]);
x=[0, 0];  % 第一个虚线的位置现在为0
y=[-1, 1];
%line(x,y,'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');

% 重新绘制第二根虚线
%t_line2 = N_precise * 0.0086 - t(drive_num+1);  % 第二条虚线的偏移后x坐标
t_line2 = 2.16;
x = [t_line2, t_line2];
y=[-1, 1];
line(x,y,'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
% text(t_line2 + 0.1, 0.5, ['$\Lambda_p = ' num2str(t_line2) '$'], ...
%     'Interpreter', 'latex', 'FontSize', 13, 'Color', 'k');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

%set(gca, 'XTick', []);  % 可根据需要调整X轴刻度
title(['$\boldmath{\omega_{19} = 1.0888,\quad \boldmath{T_s} = ' num2str(round(t_line2,2)) '}$'], ...
    'Interpreter', 'latex', 'FontSize', 20, 'Color', 'k');
ax = gca; % 获取当前轴对象
ax.XColor = 'k'; % 设置 x 轴刻度的颜色为黑色
ax.YColor = 'k'; % 设置 y 轴刻度的颜色为黑色
ax.LineWidth = 1.5; % 设置轴的边框粗细

% 如果需要单独加粗刻度的字体
ax.FontWeight = 'bold'; % 设置 x 和 y 刻度字体加粗
hold off;

 %%
axes(p(2));
% 绘制 acc_R 在偏移后的 t_offset 中的前半部分
% plot(t_offset(1:drive_num), acc_R(1:drive_num), 'Marker', 'o', 'MarkerSize', 2, ...
%     'MarkerEdgeColor', color_f1, 'MarkerFaceColor', color_f1)
% hold on;
% 
plot(t_offset(drive_num+1:timesteps), acc_R(drive_num+1:timesteps),...
    'Linestyle', '-', 'linewidth', 3, 'Color', color_f3);
hold on;
plot(t_offset(drive_num+1:timesteps), R(drive_num+1:timesteps), ...
    'Linestyle', '-.', 'linewidth', 1.5, 'Color', color_f2);
hold on;

% 设置x轴和y轴的标签
xlabel('$ {\Lambda}t $', 'interpreter', 'latex');
ylabel('\boldmath$R(t)$', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% 设置y轴范围
ylim([-0.01 0.6]);

% 设置x轴范围，偏移后的范围从0开始
%xlim([min(t_offset)-0.1 max(t_offset)+0.1-4.5]);
xlim([0 max(t_offset)+0.1-4.5]);

% 绘制第一条虚线（此时对应t_offset = 0）
x = [0, 0];  % 第一条虚线现在为0
y = [-0.01, 0.6];
%line(x, y, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');

% 绘制第二条虚线
%t_line2 = N_precise * 0.0086 - t(drive_num+1);  % 第二条虚线的偏移后x坐标
%t_line2 =2.31;
x = [t_line2, t_line2];
y = [-0.01, 0.6];
line(x, y, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
yticks([0 0.2 0.4]);
yticklabels({'0', '0.2', '0.4'});
% 添加标注 $\Lambda_p = t$，在第二条虚线附近
% text(t_line2 + 0.1, 0.1, ['$\Lambda_p = ' num2str(round(t_line2,2)) '$'], ...
%     'Interpreter', 'latex', 'FontSize', 20, 'Color','k');

% 设置图形的单位和位置
set(gcf, 'unit', 'centimeters', 'position', [8 1 12 8]);
ax = gca; % 获取当前轴对象
ax.XColor = 'k'; % 设置 x 轴刻度的颜色为黑色
ax.YColor = 'k'; % 设置 y 轴刻度的颜色为黑色
ax.LineWidth = 1.5; % 设置轴的边框粗细

% 如果需要单独加粗刻度的字体
ax.FontWeight = 'bold'; % 设置 x 和 y 刻度字体加粗

hold off;

