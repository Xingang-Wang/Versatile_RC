  function dtheta_dt = kuramoto(t, theta)
 
        dtheta_dt = omega;
        for i = 1:N
            for j = 1:N
                dtheta_dt(i) = dtheta_dt(i) + K * A(i, j) * sin(theta(j) - theta(i));
            end
        end
    end