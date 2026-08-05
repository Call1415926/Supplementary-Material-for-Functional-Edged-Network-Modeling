function output = target_f(X, X_observe_central,Omiga, U,Sigma, alpha)
%TARGET_F Objective function to be optimized
%input:X The current tensor
%     :X_observe_central The observed tensor
%     :Omiga The 0/1 projection tensor
%     :U The basis matrix to be constrained along the third mode
%     :Sigma The difference matrix for the smoothness constraint
%     :alpha The coefficient vector for the smoothness constraint
%output:The value of the objective function

part1 = 0.5 * norm(times(X - X_observe_central, Omiga))^2;
part2 = 0;
for i = 1:size(U, 2)
   part2 = part2 + 0.5 * alpha(i)* U(:,i)'* Sigma * U(:,i);
end
output = part1 + part2;
end

