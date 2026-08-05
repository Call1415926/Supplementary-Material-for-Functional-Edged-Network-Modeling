function output = mu_hat(X_observe, Omiga, seqtime)
%MU_HAT Function used to estimate mu
%input:Omiga The observation projection tensor, a 0/1 tensor with the same size as X_observe
%     :X_observe The observed tensor, with unobserved entries set to 0
%     :seqtime The full set of time points at which observations are taken
%output:The predicted values of mu(t) at all time points in seqtime, returned as a tensor whose
%       first two modes are fixed and whose third mode has the same size as seqtime, with identical values

total_obs = []; 
n = size(X_observe, 1);
m = size(X_observe, 2);
temp = double(Omiga);
N_obs = sum(temp(:)); 
for i = 1:n
   for j = 1:m
       index = find(Omiga(i, j, :) == 1);
       sub_obs = [seqtime(index); double(X_observe(i, j, index))'];
       total_obs = [total_obs, sub_obs];
   end
end
h = 1 * N_obs^(-1/5); 
Y = total_obs(2,:)';
mu = ones(size(seqtime));
for j = 1:length(seqtime)
    t = seqtime(j);
    X = ones(N_obs, 2);
    X(:, 2) = t - total_obs(1,:)';
    W = diag(ones(1,N_obs));
    for i = 1:N_obs
       W(i,i) = kappa_1D((total_obs(1,i) - t)/h); 
    end
    beta_hat = (X'* W * X) \ (X'*W*Y); 
    mu(j) = beta_hat(1);
end
mu_tensor = ones(n, m, length(mu));
for i = 1:n
    for j = 1:m
        mu_tensor(i, j, :) = mu;
    end
end
output = tensor(mu_tensor);
end

