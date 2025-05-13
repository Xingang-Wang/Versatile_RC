function u=func_get_u(theta_sin,theta_cos,rep_index,rep_num,edge,t,omega)
%Parameters
N = 50; % Number of nodes
d = 3; % Desired degree

%%
rep_all_indexs=zeros(1,d+1);
rep_all_indexs(1,1)=rep_index;
rep_all_indexs(1,2:end)=edge(edge(:,1)==rep_index,2);
 for i=1:d+1
     test_data(2*i-1)=theta_sin(rep_all_indexs(i));
     test_data(2*i)=theta_cos(rep_all_indexs(i));
 end
  test_data(2*(d+1)+1)=omega(rep_all_indexs(1,1));
  u=test_data';