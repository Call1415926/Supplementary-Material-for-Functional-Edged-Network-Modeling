function [G_new, U_new, X_new] = Retraction(X,delta,r)
%RETRACTION Symmetry-HOSVD retraction onto the FEN manifold.
%           of X + delta
%input:X The current point on the Riemannian manifold
%     :delta The increment found in the tangent space (note: X + delta is the next point)
%     :r The Tucker rank, a 1×3 vector
%output:G_new, U_new, X_new The core tensor, factor matrices, and full tensor retracted onto the
%       Riemannian manifold; U_new is a 1×3 cell array

X_trial = X + delta;

if size(X_trial,1) ~= size(X_trial,2)
    error('The first two tensor modes must have the same dimension.');
end
if r(1) ~= r(2)
    error('The first two Tucker ranks must be equal for the shared Phi factor.');
end

% Symmetry-HOSVD from Algorithm 2: the same Phi is used in modes 1 and 2.
X_mode1 = double(tenmat(X_trial,1));
X_mode2 = double(tenmat(X_trial,2));
[Phi,~,~] = svd(0.5*(X_mode1 + X_mode2),'econ');
Phi = Phi(:,1:r(1));
G_basis = nvecs(X_trial,3,r(3));

U_new = {Phi,Phi,G_basis};
G_new = ttm(X_trial,U_new,'t');
X_new = ttm(G_new,U_new);
end

