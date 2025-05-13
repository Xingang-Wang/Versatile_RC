clc;clear;
load opt_data_d_delay_resize_500_test12345.mat 
load min_rng_set_d_delay_resize_500_test12345.mat
d=2;
max_d=4;
delay=2;
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
 resSize =500;
drive_num=1000;
data_num=train_num;
testnum=size(testdata,2);
test_num=size(testdata,2);
data_len=ones(1,data_num)*train_length;
initLen = 100;
TrainLen=sum(data_len)-1;
test_Len =1000;
testLen=1000; 
outSize = 3;
nonliner_num=2;
%%
%X = zeros(nonliner_num*resSize+1,TrainLen);
X=[];
Yt = train_output_data(1:outSize,2:TrainLen+1);% run the reservoir with the data and collect X
%%
a = opt_result(3);
inSize=max_d+4;
W_in_a =opt_result(2);
Win= (2.0*rand(resSize,inSize)-1.0)*W_in_a;

%%
for i=1:train_num
    indata=train_data{i};
    [X1,W,reg]=func_get_X(resSize,opt_result,opt_rng,indata,Win,train_length);
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
%%
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

x0 = 0.01 * ones(3*N, 1);
%x0=
T=30;
tspan = [0 :0.02: 30];
dt = 0.02; % Step size for RK4

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

D=sum(A,1);
%%
rep_num=1;
% Preallocate space for solution and order parameters
num_steps = length(tspan);
X = zeros(num_steps, 3*N);
sampled_x = zeros(N,num_steps);
sampled_y = zeros(N,num_steps);
sampled_z = zeros(N,num_steps);


x1=2*rand(resSize,1)-1;
predict_time=zeros(50,N);
best=zeros(50);
%%
rng_num=5;
rep_index=6;

%%
rng(rng_num);
x0 = rand * ones(3*N, 1);
X(1, :) = x0';
y=x0;
y1=y;
combined_neighbors = cell(size(A,1),1); % 初始化合并后的邻居列表
A_squared = A^2; % 计算连接矩阵的平方
for i = 1:size(A,1)
    first_neighbors = find(A(i,:) > 0); % 第一邻居
    neighbors_me{i}=[i,first_neighbors];
end
for idx = 1:num_steps-1
    t = tspan(idx);
    %y = X(idx, :)';
    
    k1 = dt * coupled_lorenz_multi_adj(t, y, sigmas, betas, rhos, kappa, N, A);
    k2 = dt * coupled_lorenz_multi_adj(t + dt/2, y + k1/2, sigmas, betas, rhos, kappa, N, A);
    k3 = dt * coupled_lorenz_multi_adj(t + dt/2, y + k2/2, sigmas, betas, rhos, kappa, N, A);
    k4 = dt * coupled_lorenz_multi_adj(t + dt, y + k3, sigmas, betas, rhos, kappa, N, A);
    
    y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
    X(idx+1, :) = y';
    
    for i = 1:N
     sampled_x(i,:) = X(:,3*i-2)';
     sampled_y(i,:) = X(:,3*i-1)';
     sampled_z(i,:) = X(:,3*i)';
     sampled_x(i,:) = 2*((sampled_x(i,:)-minx)./(maxx-minx))-1;
     sampled_y(i,:)= 2*((sampled_y(i,:)-miny)./(maxy-miny))-1;
    sampled_z(i,:)= 2*((sampled_z(i,:)-minz)./(maxz-minz))-1;
    end
    % sampled_x = 2*((sampled_x-minx)./(maxx-minx))-1;
    % sampled_y= 2*((sampled_y-miny)./(maxy-miny))-1;
    % sampled_z= 2*((sampled_z-minz)./(maxz-minz))-1;
    % 

    if idx>delay
    test_data=zeros(max_d+4,1);
    num=size(neighbors_me{rep_index},2);
    for j=1:num-1
        Labels=neighbors_me{rep_index}';
    test_data(j)=sampled_x(Labels(j+1),idx)-sampled_x(Labels(1),idx);
    end
    for k=num:max_d
        test_data(k)=sampled_x(Labels(2),idx+(k-(max_d+1))+delay);
    end
    test_data(max_d+1)=sampled_x(Labels(1),idx);
    test_data(max_d+2)=sampled_y(Labels(1),idx);
    test_data(max_d+3)=sampled_z(Labels(1),idx);
    test_data(max_d+4,:)=c3s(Labels(1))/50;
      u=test_data;
     x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
     RC_p = Wout*[1;x1;x1.^2;];
     RC_x=RC_p(1);
     RC_y=RC_p(2);
     RC_z=RC_p(3);
     
     RC_rep_x = ((RC_x+1)/2)*(maxx-minx)+minx;
     RC_rep_y = ((RC_y+1)/2)*(maxy-miny)+miny;
     RC_rep_z = ((RC_z+1)/2)*(maxz-minz)+minz;
     
    if idx>drive_num
        y(3*rep_index-2)= RC_rep_x;
        y(3*rep_index-1)= RC_rep_y;
        y(3*rep_index)= RC_rep_z;
    end
    y1=y;
    y1(3*rep_index-2)= RC_rep_x;
    y1(3*rep_index-1)= RC_rep_y;
    y1(3*rep_index)= RC_rep_z;
    end
    
    pre_x2(idx)=y(4);
    pre_xrep(idx)=y(3*rep_index-2);
   %%
    % Calculate average Euclidean distance for order parameter
    total_distance = 0;
    count = 0;
    for i = 1:N
        for j = i+1:N
            xi = y1(3*i-2:3*i);
            xj = y1(3*j-2:3*j);
            distance = norm(xi - xj);
            total_distance = total_distance + distance;
            count = count + 1;
        end
    end
    R(idx+1) = total_distance / count; % Store the average distance
end
[acc_R,acc_X]=func_acc_R_loz(rng_num);
%%
color_f1= addcolorplus(16);
color_f3= addcolorplus(3);
color_f2= addcolorplus(185);
%color_f2= addcolorplus(73);
%color_f4= addcolorplus(215);
%%
true_drive_num=870;
drive_num=50;

pre_xrep = pre_xrep(true_drive_num+1:end);
acc_X = acc_X(true_drive_num+1:end,:);
R=R(true_drive_num+1:end);
acc_R=acc_R(true_drive_num+1:end,:);
%%
epsilon_percent=0.1;
tau=20;
N_precise = func_calculate_precise_prediction_periods(acc_R, R, epsilon_percent, tau);
%N_precise=N_precise(true_drive_num+1:end);
%N_precise=N_precise-true_drive_num;
le=1.09;
plot_drive=50;
t = (0.02:0.02:length(R)*0.02)*le;
t_offset = t- t(drive_num+1);
%%
TightPlot.ColumeNumber = 2;     % 子图行数
TightPlot.RowNumber = 1;    % 子图列数
TightPlot.GapW = 0.2;  % 子图之间的左右间距
TightPlot.GapH = 0.0;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.2;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.13;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.16;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.015;  % 子图与图片右方的间距


p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]); 
%%
t_line2 =3.05;
figure
axes(p(1));
t_offset = t - t(drive_num+1);  % 偏移 t，使虚线位置为 0

% 重新绘制图像
% plot(t_offset(1:drive_num),acc_X(1:drive_num,3*rep_index-2),'Marker','o', 'MarkerSize', 2, ...
%    'MarkerEdgeColor', color_f1,'MarkerFaceColor', color_f1,'Linestyle', 'none');

%hold on;
plot(t_offset(drive_num+1:end-1),acc_X(drive_num+1:length(t)-1,3*rep_index-2),...
   'Linestyle','-','linewidth',3,'Color',color_f3);
hold on;
plot(t_offset(drive_num+1:end-1),pre_xrep(drive_num+1:length(t)-1),...
   'Linestyle','-.','linewidth',1.5,'Color',color_f2);

ylim([-19 20]);
xlabel('$ \Lambda $','interpreter','latex');
ylabel('\boldmath$ x(t)  $','interpreter','latex');
hold on;
title(['$\mathbf{\rho_{6} = 33.60, T_s = ' num2str(round(t_line2,2)) ', k=1}$'], ...
    'Interpreter', 'latex', 'FontSize', 20, 'Color', 'k');

ax = gca; % 获取当前轴对象
ax.XColor = 'k'; % 设置 x 轴刻度的颜色为黑色
ax.YColor = 'k'; % 设置 y 轴刻度的颜色为黑色
ax.LineWidth = 1.5; % 设置轴的边框粗细

% 如果需要单独加粗刻度的字体
ax.FontWeight = 'bold'; % 设置 x 和 y 刻度字体加粗
% 调整 xlim，以虚线位置为0，范围重新设置为 [0, 原始t的最大值 - 虚线位置对应的t值]
%xlim([min(t_offset)-0.1 max(t_offset)+0.1-5]);
xlim([0 max(t_offset)+0.1-5]);
% x=[0, 0];  % 第一个虚线的位置现在为0
% y=[-19, 20];
% line(x,y,'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
%t_line2 = N_precise *0.02*le - t(drive_num+1);  % 第二条虚线的偏移后x坐标

% 重新绘制第二根虚线
x=[t_line2, t_line2];
y=[-19, 20];
line(x,y,'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');
% text(t_line2 + 0.1, 0.5, ['$\Lambda_p = ' num2str(t_line2) '$'], ...
%     'Interpreter', 'latex', 'FontSize', 13, 'Color', 'k');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
yticks([-10 5 20]);
yticklabels({'-10', '5', '20'});

set(gca, 'XTick', []);  % 可根据需要调整X轴刻度
%hold off

axes(p(2));
% 绘制 acc_R 在偏移后的 t_offset 中的前半部分
plot(t_offset(1:drive_num), acc_R(1:drive_num), 'Marker', 'o', 'MarkerSize', 2, ...
    'MarkerEdgeColor', color_f1, 'MarkerFaceColor', color_f1,'Linestyle', 'none');
hold on;

plot(t_offset(drive_num+1:end), acc_R(drive_num+1:end),...
    'Linestyle', '-', 'linewidth', 3, 'Color', color_f3);
hold on;
plot(t_offset(drive_num+1:end), R(drive_num+1:end), ...
    'Linestyle', '-.', 'linewidth', 1.5, 'Color', color_f2);
hold on;

% 设置x轴和y轴的标签
xlabel('$ \Lambda $t', 'interpreter', 'latex');
ax = gca; % 获取当前轴对象
ax.XColor = 'k'; % 设置 x 轴刻度的颜色为黑色
ax.YColor = 'k'; % 设置 y 轴刻度的颜色为黑色
ax.LineWidth = 1.5; % 设置轴的边框粗细

% 如果需要单独加粗刻度的字体
ax.FontWeight = 'bold'; % 设置 x 和 y 刻度字体加粗
ylabel('\boldmath$\delta e(t)$', 'Interpreter', 'latex');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);

% 设置y轴范围
ylim([10 31]);

% 设置x轴范围，偏移后的范围从0开始
xlim([0 max(t_offset)+0.1-5]);

% 绘制第一条虚线（此时对应t_offset = 0）
% x = [0, 0];  % 第一条虚线现在为0
% y = [10 31];
% line(x, y, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');

% 绘制第二条虚线
%t_line2 = N_precise * 0.02*le- t(drive_num+1);  % 第二条虚线的偏移后x坐标
x = [t_line2, t_line2];
y = [10, 31];
line(x, y, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k');

% 添加标注 $\Lambda_p = t$，在第二条虚线附近
% text(t_line2 + 0.1, 10, ['$\Lambda_p = ' num2str(t_line2) '$'], ...
%     'Interpreter', 'latex', 'FontSize', 20, 'Color','k');

% 设置图形的单位和位置
set(gcf, 'unit', 'centimeters', 'position', [8 1 12 8]);
hold off;
%%

 function dx = coupled_lorenz_multi_adj(t, x, sigmas, betas, rhos, kappa, N, A)
        dx = zeros(3 * N, 1);
        for i = 1:N
            xi = x(3 * i - 2:3 * i);
            
            % 洛伦兹系统方程
            dx1dt = sigmas(i) * (xi(2) - xi(1));
            dy1dt = xi(1) * (rhos(i) - xi(3)) - xi(2);
            dz1dt = xi(1) * xi(2) - betas(i) * xi(3);
            
            % 耦合项
            coupling = 0;
            for j = 1:N
                if A(i, j) > 0
                    xj = x(3 * j - 2);
                    coupling = coupling + kappa * (xj - xi(1));
                end
            end
            dx(3 * i - 2:3 * i) = [dx1dt + coupling; dy1dt; dz1dt];
        end
    end
