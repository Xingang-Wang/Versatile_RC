function N_precise = func_calculate_precise_prediction_periods(y_true, y_pred, epsilon_percent, tau)
    % 初始化变量
    N_precise = 0;
    N_exceed = 0;
    
    % 获取数据长度
    n = length(y_true);
    
    % 遍历误差序列
    for i = 1:n
        % 计算当前时间点的误差阈值
        epsilon_i = epsilon_percent * y_true(i);
        
        % 计算当前时间点的误差
        error = abs(y_true(i) - y_pred(i));
        
        if error <= epsilon_i
            % 误差小于等于当前误差阈值，重置N_exceed并增加N_precise
            N_exceed = 0;
            N_precise = N_precise + 1;
        else
            % 误差大于当前误差阈值，增加N_exceed
            N_exceed = N_exceed + 1;
            
            if N_exceed > tau
                % 当连续超出误差阈值的时间超过持续时间阈值时，停止统计
                break;
            end
        end
    end
end