function [G_next,U_next,X_next,score_next, score_list] = conjugate_dradient(G_old, U_old, X_old,eta_old,...
    G_new, U_new, X_new,score_new,X_observe_central, Omiga, Sigma, Alpha,alpha,max_m,r, tol,maxiter)
%CONJUGATE_GRADIENT Method for solving using the conjugate gradient method
%input:G_old Initial Tucker core tensor
%     :U_old Initial Tucker factor matrices, a 1×3 cell
%     :X_old Initial predicted tensor, the product of G_old and U_old
%     :eta_old Initial descent direction
%     :G_new Tucker core tensor after one iteration
%     :U_new Tucker factor matrices after one iteration, a 1×3 cell
%     :X_new Predicted tensor after one iteration, the product of G_new and U_new
%     :score_new Score after one iteration
%     :X_observe_central Centered observed tensor
%     :Omiga Mask tensor with 0/1 values
%     :Sigma Difference matrix for smoothness constraint
%     :Alpha Coefficient matrix for smoothness constraint
%     :alpha Coefficient vector for smoothness constraint
%     :r Given Tucker decomposition rank, a 1×3 vector
%     :max_m Maximum backtracking steps m in conjugate gradient
%     :tol Tolerance threshold for convergence

%output:G_next Converged Tucker core tensor
%      :U_next Converged Tucker factor matrices, a 1×3 cell
%      :X_next Converged predicted tensor, the product of G_next and U_next
%      :score_next Converged final score
%      :score_list List recording the score at each iteration
score_list = score_new;
count = 1;
while true
   grad_Euc_new = grad_Euc(G_new, U_new, X_new,X_observe_central, Omiga, Sigma, Alpha);
   epsilon_new = Euc_to_Rie(grad_Euc_new, G_new, U_new);
   beta_new = beta(G_old, U_old, X_old,G_new, U_new, X_new,X_observe_central, Omiga, Sigma, Alpha);
   vec_trans = vector_transport(eta_old, G_new, U_new);
   eta_new = -epsilon_new + beta_new * vec_trans;
   alphA = find_alpha(X_new, eta_new,Omiga,X_observe_central);
   m = find_m(X_new, alphA, U_new,eta_new, epsilon_new,r, X_observe_central,Omiga,Sigma, alpha, max_m);
   [G_next, U_next, X_next] = Retraction(X_new, 2^(-m) * alphA * eta_new, r);
   score_next = target_f(X_next, X_observe_central,Omiga, U_next{3},Sigma, alpha); 
   score_list = [score_list, score_next];
   if abs(score_next - score_new) <= tol 
      break; 
   end
   G_old = G_new;
   U_old = U_new;
   X_old = X_new;
   eta_old = eta_new;
   
   G_new = G_next;
   U_new = U_next;
   X_new = X_next;
   score_new = score_next;
   
   count = count + 1;
   if count > maxiter
      break; 
   end
end

end

