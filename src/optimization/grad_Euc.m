function output = grad_Euc(G, U, X_pre,X_observe_central, Omiga, Sigma, Alpha)
%GRAD_EU Used to compute, at X_pre, the Euclidean gradient of the objective function f
%        after adding a smoothness constraint on the third mode
%input:G The Tucker core tensor
%     :U The Tucker factor matrices, a 1×3 cell
%     :X_pre The tensor obtained by multiplying G and U in the Tucker decomposition
%     :X_observe_central The observed tensor (centered)
%     :Omiga The projection mask tensor for partial observations
%     :Sigma The difference matrix for the smoothness constraint
%     :Alpha The diagonal coefficient matrix for the smoothness constraint
%output:The Euclidean gradient of f at X_pre after adding the smoothness constraint on the third mode,
%       a tensor with the same size as the observed tensor

grad_1 = times(X_pre - X_observe_central, Omiga);
G_mat = ttm(G, {U{1}, U{2}}, [1,2]);
G_mat = double(tenmat(G_mat, 3));
S = G_mat'/(G_mat*G_mat');
grad_2 = Sigma * U{3} * Alpha *S';
temp = tenmat(X_observe_central, 3);
grad_2 = tenmat(grad_2, temp.rdims, temp.cdims, temp.tsize); 
grad_2 = tensor(grad_2);
output = grad_1 + grad_2;
end

