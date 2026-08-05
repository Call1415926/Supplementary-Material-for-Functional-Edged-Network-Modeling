function output = beta(G_old, U_old, X_old,G_new, U_new, X_new,X_observe_central, Omiga, Sigma, Alpha)

grad_Euc_old = grad_Euc(G_old, U_old, X_old,X_observe_central, Omiga, Sigma, Alpha);
grad_Euc_new = grad_Euc(G_new, U_new, X_new,X_observe_central, Omiga, Sigma, Alpha);
grad_Rie_old = Euc_to_Rie(grad_Euc_old, G_old, U_old);
grad_Rie_new = Euc_to_Rie(grad_Euc_new, G_new, U_new);
trans_vec = vector_transport(grad_Rie_old,  G_new, U_new);
output = max(0, innerP(grad_Rie_new, grad_Rie_new - trans_vec)/norm(grad_Rie_old)^2);


end

