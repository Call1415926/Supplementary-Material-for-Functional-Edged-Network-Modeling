function [X_output, G_output, U_output] = adjust(X_in, G_in, U_in, X_true, G_true, U_true,t_start, t_end, seq)
%ADJUST Orthogonally align estimated factors for reporting only.
% The transformation is orthogonal and therefore does not change X_in.
% X_true and G_true are retained in the interface for backward compatibility.
seqtime = t_start:seq:t_end;
constant = sqrt(length(seqtime)/2); 
U_true{3} = U_true{3}/constant;

% A single Procrustes rotation is used for the shared node factor.
[L,~,R] = svd(U_in{1}'*U_true{1},'econ');
Q_Phi = L*R';

% An independent orthogonal rotation is permitted in the functional mode.
[L,~,R] = svd(U_in{3}'*U_true{3},'econ');
Q_G = L*R';

O = {Q_Phi,Q_Phi,Q_G};
U_output = {U_in{1}*Q_Phi,U_in{2}*Q_Phi,U_in{3}*Q_G};
G_output = ttm(G_in, O, 't');
X_output = X_in;
end

