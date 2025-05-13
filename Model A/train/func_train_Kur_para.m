function [rmse,rmse_dynamic]  = func_train_Kur_para(hyperpara_set,rng_num,train_data,train_output_data,testdata,test_output_data,train_length,train_num)
%%
d=3;
max_d=3;
rng(rng_num);
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
W_in_a = hyperpara_set(2);
Win= (2.0*rand(resSize,inSize)-1.0)*W_in_a;


%%
for i=1:train_num
    indata=train_data{i};
    [X1,W,reg]=func_get_X(inSize,hyperpara_set,rng_num,indata,Win,train_length);
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
Y1= zeros(outSize,testLen);
a = hyperpara_set(3);
rmse_dynamic=0;  
rng(rng_num);

for i=1:testnum
    Testdata=testdata{i};
    u=Testdata(1:2*(d+1)+1,1);
    Y1(:,1)=u(1:2);
    x1=2*rand(resSize,1)-1;
for t = 1:test_Len-1 
    x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
    y = Wout*[1;x1;x1.^2;];
    Y1(:,t+1) = y;
     if t<drive_num
     u(1)=Testdata(1,1+t);
     u(2)=Testdata(2,1+t);
    else
    u(1:2) = y;   
     end  
    u(3:end)=Testdata(3:end,1+t);
    
end
rmse_dynamic=rmse_dynamic+mean(abs(Y1(1,drive_num+1:testLen)-...
    Testdata(1,drive_num+1:testLen)))+...
    mean(abs(Y1(2,drive_num+1:testLen)-...
    Testdata(2,drive_num+1:testLen)));
end
rmse_dynamic=rmse_dynamic/testnum;
rmse=rmse_dynamic;
if isnan(rmse) || rmse>10
    rmse=10;
end
