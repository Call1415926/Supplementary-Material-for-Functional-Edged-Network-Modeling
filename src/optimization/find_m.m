function output = find_m(X_new, alphA, U_new,eta_new, epsilon_new,r,...
    X_observe_central,Omiga,Sigma, alpha, max_m)
%FIND_M Find the optimal backtracking length m
%input:X_new The current predicted tensor
%     :X_observe_central The observed tensor
%     :Omiga The 0/1 projection mask tensor
%     :U_new The current Tucker factor matrices, a 1×3 cell
%     :Sigma The difference matrix for the smoothness constraint
%     :alpha The coefficient vector for the smoothness constraint
%     :alphA The currently computed optimal step size
%     :eta_new The currently computed descent direction
%     :epsilon_new The Riemannian gradient at the current point X_new
%     :r The Tucker decomposition rank, a 1×3 vector
%     :max_m The maximum allowed m for backtracking; if no suitable m is found before max_m, use max_m
%output:The optimal backtracking step length m

for m = 0:max_m
    score_1 = target_f(X_new, X_observe_central,Omiga, U_new{3},Sigma, alpha);%计算第一个得分
    [~, U_other, X_other] = Retraction(X_new,2^(-m)*alphA*eta_new,r);
    score_2 = target_f(X_other, X_observe_central,Omiga, U_other{3},Sigma, alpha);%计算第二个得分，也即下降后的得分
    my_inner = innerP(epsilon_new, 2^(-m) * alphA * eta_new);
    if score_1 - score_2 >= - 10^(-4) * my_inner
       break; 
    end
end
output = m;

end

