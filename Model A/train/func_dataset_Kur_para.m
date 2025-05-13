function [train_data,train_output_data,testdata,test_output_data,train_length,train_num]=func_dataset_Kur_para()
%clc;clear
train_length=4000;
test_length=2000;
[theta_sin,theta_cos,omega,A,d]=evo_kur();
max_d=max(d);
%save kur_N_50.mat  omega A d
interval=10;
delay=0;
sampled_theta_sin = [];
sampled_theta_cos = [];
for i = 1:size(theta_sin, 1)
    sampled_row = downsample(theta_sin(i,:), interval);
    sampled_theta_sin = [sampled_theta_sin; sampled_row];
    sampled_row = downsample(theta_cos(i,:), interval);
    sampled_theta_cos = [sampled_theta_cos; sampled_row];
end
% 计算邻居和第二邻居
% 计算邻居和第二邻居
neighbors_me = cell(size(A,1),1); % 初始化合并后的邻居列表
A_squared = A^2; % 计算连接矩阵的平方
for i = 1:size(A,1)
    first_neighbors = find(A(i,:) > 0); % 第一邻居
    neighbors_me{i}=[i,first_neighbors];
end
%%

rng(1);
index_num=10;
train_test_index=randperm(50,40);
train_index=train_test_index(1:index_num);
%save train_index_10_rng_5.mat train_index
train_num=size(train_index,2);
Train_data=cell(1,train_num);
%%
transient=1001;
train_output_data=zeros(2,train_num*train_length);
for j=1:train_num
    train_data=zeros(2*(max_d+1)+1,train_length);
   
    for i=1:size(neighbors_me{train_index(j)},2)
        Labels=neighbors_me{train_index(j)}';
        train_data(2*i-1,:)=sampled_theta_sin(Labels(i),transient+delay:...
            transient+train_length-1+delay);

        train_data(2*i,:)=sampled_theta_cos(Labels(i),transient+delay:...
            transient+train_length-1+delay);
    end
       train_output_data(1,(j-1)*train_length+1:j*train_length)=sampled_theta_sin(Labels(1),transient+delay:...
            transient+train_length-1+delay);
        train_output_data(2,(j-1)*train_length+1:j*train_length)=sampled_theta_cos(Labels(1),transient+delay:...
            transient+train_length-1+delay);
        train_data(2*(max_d+1)+1,:)=omega(Labels(1))*ones(1,train_length);
        Train_data{j}=train_data;
end
%%'
transient=5001;
test_index=train_index;
test_num=size(test_index,2);
Testdata=cell(1,test_num);

for i=1:test_num
    testdata=zeros(size(neighbors_me{test_index(i)},2)+2,test_length);
    test_output_data=zeros(test_num*2,test_length);
    for j=1:1:size(neighbors_me{test_index(i)},2)
        Labels=neighbors_me{test_index(i)}';
    testdata(2*j-1,:)=sampled_theta_sin(Labels(j),transient+delay:transient+test_length-1+delay);
    testdata(2*j,:)=sampled_theta_cos(Labels(j),transient+delay:transient+test_length-1+delay);

    test_output_data(2*i-1,:)=sampled_theta_sin(Labels(1),transient+train_length:transient+train_length+test_length-1);
    test_output_data(2*i+1,:)=sampled_theta_cos(Labels(1),transient+train_length:transient+train_length+test_length-1);
    end
     testdata(2*(max_d+1)+1,:)=omega(Labels(1))*ones(1,test_length);
      Testdata{i}=testdata;
end
train_data=Train_data;
testdata=Testdata;
% figure
% for i=2:4
%     plot(testdata(2*i,1:1000),'.-');
%     hold on;
% end