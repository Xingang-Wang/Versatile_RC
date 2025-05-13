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
