function output = find_alpha(X_new, eta_new,Omiga,X_observe_central)
%FIND_ALPHA Used to compute the optimal alphA in the conjugate gradient algorithm, i.e., the step size
%input:X_new The current estimated tensor
%     :eta_new The current descent direction
%     :Omiga The mask 0/1 tensor
%     :X_oberve_central The observed tensor
%output:The optimal step size alphA (a scalar). Note: distinguish this from alpha in the smoothness constraint;
%       here the last letter is capitalized on purpose.

inner_1 = innerP(times(Omiga, eta_new), times(Omiga, X_observe_central-X_new));
inner_2 = norm(times(Omiga, eta_new))^2;
if inner_2 <= eps
    output = 0;
else
    output = inner_1/inner_2;
end
if ~isfinite(output) || output <= 0
    % The closed-form value uses the data-fit term only. If the smoothness
    % term makes it nonpositive, start Armijo backtracking from one.
    output = 1;
end
end

