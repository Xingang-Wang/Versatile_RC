function Testdata=func_testdata_chua(minx,maxx,miny,maxy,minz,maxz)
test_length=3000;
train_length=3000;
[X,A,N,c3s]=func_chua_evolve(minx,maxx,miny,maxy,minz,maxz);
d=sum(A,1);
max_d=max(d);
%delay=max(d)-min(d);
delay=2;
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
test_index=[1,2,3,4,5];
test_num=size(test_index,2);
Testdata=cell(1,test_num);
transient=1000;
for i=1:test_num
    testdata=zeros(3*(max_d+1+delay)+1,test_length);
    test_output_data=zeros(test_num*2,test_length);
    for j=1:1:size(neighbors_me{test_index(i)},2)-1
        Labels=neighbors_me{test_index(i)}';
    testdata(3*j-2,:)=sampled_x(Labels(j+1),transient+delay:transient+test_length-1+delay);
    testdata(3*j-1,:)=sampled_y(Labels(j+1),transient+delay:transient+test_length-1+delay);
    testdata(3*j,:)=sampled_z(Labels(j+1),transient+delay:transient+test_length-1+delay);
    test_output_data(2*i-1,:)=sampled_x(Labels(1),transient+train_length:transient+train_length+test_length-1);
    test_output_data(2*i+1,:)=sampled_y(Labels(1),transient+train_length:transient+train_length+test_length-1);
    end
    for k=size(neighbors_me{test_index(i)},2):size(neighbors_me{test_index(i)},2)+delay
        testdata(3*k-2,:)=sampled_x(Labels(1),transient+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay:...
            transient+test_length-1+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay);
        testdata(3*k-1,:)=sampled_y(Labels(1),transient+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay:...
            transient+test_length-1+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay);
        testdata(3*k,:)=sampled_z(Labels(1),transient+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay:...
            transient+test_length-1+(k-(size(neighbors_me{test_index(i)},2)+delay))+delay);
    end
      testdata(3*(max_d+1+delay)+1,:)=c3s(Labels(1))*ones(1,test_length)/50;
      %estdata(3*(max_d+1)+2,:)=(size(neighbors_me{test_index(i)},2)-1)/max_d*ones(1,test_length);
      Testdata{i}=testdata;
end