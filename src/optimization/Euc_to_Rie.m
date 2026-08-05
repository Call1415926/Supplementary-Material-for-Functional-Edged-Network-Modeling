function output = Euc_to_Rie(A, G, U)
%EUC_TO_RIE Project A onto the shared-Phi FEN tangent space.
%input:A The target tensor to be projected
%     :G The Tucker core at point X_pre
%     :U A 1×3 cell array of Tucker factor matrices at point X_pre
%     :X_pre The tangent point on the Riemannian manifold for projection
%output:The projection result, which is a tensor


if size(U{1},1) ~= size(U{2},1) || size(U{1},2) ~= size(U{2},2)
    error('The first two factors must have matching dimensions.');
end
if norm(U{1}-U{2},'fro') > 1e-8
    error('The FEN manifold requires U{1} and U{2} to be the same Phi.');
end

Phi = U{1};
G_basis = U{3};

% Candidate node-factor tangent variations from modes 1 and 2.
temp1 = eye(size(Phi,1)) - Phi*Phi';
temp2 = double(tenmat(ttm(A, {U{2}, U{3}}, [2, 3], 't'), 1));
temp = double(tenmat(G, 1));
temp3 = temp' * pinv(temp*temp');
W1 = temp1*temp2*temp3;

temp2 = double(tenmat(ttm(A, {U{1}, U{3}}, [1, 3], 't'), 2));
temp = double(tenmat(G, 2));
temp3 = temp' * pinv(temp*temp');
W2 = temp1*temp2*temp3;

% One shared variation is used in both node modes, preserving U1 = U2.
W_Phi = 0.5*(W1 + W2);

temp1 = eye(size(G_basis,1)) - G_basis*G_basis';
temp2 = double(tenmat(ttm(A, {U{1}, U{2}}, [1, 2], 't'), 3));
temp = double(tenmat(G, 3));
temp3 = temp' * pinv(temp*temp');
W_G = temp1*temp2*temp3;

output = ttm(A,{Phi*Phi',Phi*Phi',G_basis*G_basis'}) + ...
    ttm(G,{W_Phi,Phi,G_basis}) + ...
    ttm(G,{Phi,W_Phi,G_basis}) + ...
    ttm(G,{Phi,Phi,W_G});

end

