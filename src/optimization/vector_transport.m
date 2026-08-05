function output = vector_transport(epsilon_old, G_new, U_new)
%VECTOR_TRANSPORT y is used to compute the transported vector in the conjugate gradient iteration,
%                 transported from epsilon_old
%input:elpsion_old The descent tensor from the previous iteration, or any tensor to be transported
%     :G_new The current Tucker core tensor
%     :U_new The current Tucker factor matrices, a 1×3 cell
%output:The transported tensor, which can be used for the next conjugate update

% Projection transport is the tangent projection at the new point. Calling
% the same shared-Phi projector ensures the transported direction remains in
% the FEN tangent space.
output = Euc_to_Rie(epsilon_old,G_new,U_new);
end

