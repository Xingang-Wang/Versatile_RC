
function [A,omega]=func_connectivity_A(N,gamma,delta,d)
rng(2);
%Parameters
% N = 50; % Number of nodes
% d = 3; % Desired degree
% gamma = 5; % Assortativity parameter
% delta=0.8;
% Initialize node frequencies
omega = rand(N, 1) * pi - pi/2;

% Initialize adjacency matrix
A = zeros(N,N);

% Start adding links
for i = 1:N
    % Check if node i needs additional links
    while sum(A(i, :)) < d
        % Randomly pick another node j
        j = randi(N);
        % Ensure j is not already connected to i and requires additional links
        while j == i || A(i, j) == 1 || sum(A(j, :)) >= d
            j = randi(N);
        end
        % Calculate the linking probability
        p_ij = delta^gamma/ (delta^gamma + abs(omega(i) - omega(j))^gamma);
        % Link nodes i and j with probability p_ij
        if rand() < p_ij
            A(i, j) = 1;
            A(j, i) = 1; % Assuming undirected network
        end
    end
end


