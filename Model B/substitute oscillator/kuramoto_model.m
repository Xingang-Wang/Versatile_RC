function dtheta_dt = kuramoto_model(theta, omega, A, K)
    N = length(theta);
    dtheta_dt = omega + K*sum(A.*sin(ones(N, 1) * theta'-theta * ones(1, N)),2);
end