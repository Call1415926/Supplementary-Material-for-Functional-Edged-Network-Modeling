function [G,U,X_pre] = Initialization_tucker(X, r)
%INITIALIZATION_TUCKER Symmetry-preserving Tucker initialization.
%input:X The input observed tensor, i.e., the centered discrete observations
%     :r A 1×3 vector specifying the Tucker rank in each of the three modes
%output:G The initial Tucker core tensor
%      :U A 1×3 cell array, where each element is an initial Tucker factor matrix
%      :X_pre The tensor obtained by multiplying G and U

if size(X,1) ~= size(X,2)
    error('The first two tensor modes must have the same dimension.');
end
if r(1) ~= r(2)
    error('The first two Tucker ranks must be equal for the shared Phi factor.');
end

% Algorithm 2 in the paper: estimate one common node factor Phi from the
% two node-mode matricizations and use it in both modes.
X_mode1 = double(tenmat(X,1));
X_mode2 = double(tenmat(X,2));
[Phi,~,~] = svd(0.5*(X_mode1 + X_mode2),'econ');
Phi = Phi(:,1:r(1));

% Dominant left singular vectors in the functional mode.
G_basis = nvecs(X,3,r(3));
U = {Phi,Phi,G_basis};
G = ttm(X, U, 't');
X_pre = ttm(G, U);
end

