%function [train_data,train_output_data,testdata,test_output_data,train_length,train_num,minx,maxx,miny,maxy,minz,maxz]=func_dataset_chua_classic_coupled()
clc;clear
train_length=3000;
test_length=3000;
%[theta_sin,theta_cos,omega,A,d]=evo_kur_d();
[X,A,N,c3s,minx,maxx,miny,maxy,minz,maxz]=func_chua_classic_coupled();
d=sum(A,1);
max_d=max(d);
delay=max(d)-min(d);
%save kur_dataset_A_rng2.mat data_theta omega A d
interval=20;
sampled_x= [];
sampled_y= [];
sampled_z= [];
for i = 1:N 
     sampled_row = downsample(X(:,3*i-2) ,interval);
     sampled_x = [sampled_x;sampled_row'];
    sampled_row = downsample(X(:,3*i-1), interval);
    sampled_y = [sampled_y;sampled_row'];
    sampled_row = downsample(X(:,3*i), interval);
    sampled_z = [sampled_z;sampled_row'];
end


% 计算邻居和第二邻居
neighbors_me = cell(size(A,1),1); % 初始化合并后的邻居列表
A_squared = A^2; % 计算连接矩阵的平方
for i = 1:size(A,1)
    first_neighbors = find(A(i,:) > 0); % 第一邻居
    neighbors_me{i}=[i,first_neighbors];
end
%%
rng(1);
index_num=4;
train_test_index=randperm(N,index_num);
%train_index=train_test_index(1:index_num/2);
train_index=[1,2,3];
train_num=size(train_index,2);
Train_data=cell(1,train_num);
%%
transient=2001;
train_output_data=zeros(3,train_num*train_length);
for j=1:train_num
    train_data=zeros(max_d+4,train_length);
    num=size(neighbors_me{train_index(j)},2);
    for i=1:size(neighbors_me{train_index(j)},2)-1
        Labels=neighbors_me{train_index(j)}';
        train_data(i,:)=sampled_x(Labels(i+1),transient+delay:...
            transient+train_length-1+delay)-sampled_x(Labels(1),transient+delay:...
            transient+train_length-1+delay);
    end
    for k=num:max_d
        train_data(k,:)=sampled_x(Labels(2),transient+(k-(max_d+1))+delay:...
            transient+train_length-1+(k-(max_d+1))+delay)-...
            sampled_x(Labels(1),transient+(k-(max_d+1))+delay:...
            transient+train_length-1+(k-(max_d+1))+delay);
    end
   train_data(max_d+1,:)=sampled_x(Labels(1),transient+delay:...
            transient+train_length-1+delay);
   train_data(max_d+2,:)=sampled_y(Labels(1),transient+delay:...
            transient+train_length-1+delay);
   train_data(max_d+3,:)=sampled_z(Labels(1),transient+delay:...
            transient+train_length-1+delay);
   
       train_output_data(1,(j-1)*train_length+1:j*train_length)=sampled_x(Labels(1),transient+delay:...
            transient+train_length-1+delay);
        train_output_data(2,(j-1)*train_length+1:j*train_length)=sampled_y(Labels(1),transient+delay:...
            transient+train_length-1+delay);
        train_output_data(3,(j-1)*train_length+1:j*train_length)=sampled_z(Labels(1),transient+delay:...
            transient+train_length-1+delay);
        train_data(max_d+4,:)=c3s(Labels(1))*ones(1,train_length)/50;
        %train_data(3*(max_d+1)+2,:)=(size(neighbors_me{train_index(j)},2)-1)/max_d*ones(1,train_length);
        Train_data{j}=train_data;
      
end
%%
%test_index=train_test_index(index_num/2+1:end);
%test_index=train_index;
test_index=[1,2,3,4,5];
test_num=size(test_index,2);
Testdata=cell(1,test_num);
transient=5001;
for i=1:test_num
    testdata=zeros(max_d+4,test_length);
    test_output_data=zeros(test_num*2,test_length);
    num=size(neighbors_me{test_index(i)},2);
    for j=1:size(neighbors_me{test_index(i)},2)-1
        Labels=neighbors_me{test_index(i)}';
    testdata(j,:)=sampled_x(Labels(j+1),transient+delay:transient+test_length-1+delay)-sampled_x(Labels(1),transient+delay:...
            transient+test_length-1+delay);
    test_output_data(2*i-1,:)=sampled_x(Labels(1),transient+train_length:transient+train_length+test_length-1);
    test_output_data(2*i+1,:)=sampled_y(Labels(1),transient+train_length:transient+train_length+test_length-1);
    end
    for k=num:max_d
        testdata(k,:)=sampled_x(Labels(2),transient+(k-(max_d+1))+delay:...
            transient+train_length-1+(k-(max_d+1))+delay)-...
            sampled_x(Labels(1),transient+(k-(max_d+1))+delay:...
            transient+train_length-1+(k-(max_d+1))+delay);
    end
   testdata(max_d+1,:)=sampled_x(Labels(1),transient+delay:...
            transient+test_length-1+delay);
   testdata(max_d+2,:)=sampled_y(Labels(1),transient+delay:...
            transient+test_length-1+delay);
   testdata(max_d+3,:)=sampled_z(Labels(1),transient+delay:...
            transient+test_length-1+delay);
      testdata(max_d+4,:)=c3s(Labels(1))*ones(1,test_length)/50;
      %estdata(3*(max_d+1)+2,:)=(size(neighbors_me{test_index(i)},2)-1)/max_d*ones(1,test_length);
      Testdata{i}=testdata;
end
train_data=Train_data;
testdata=Testdata;

