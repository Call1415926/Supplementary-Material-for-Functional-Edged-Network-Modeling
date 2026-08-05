function [MSE,R2,Rec_esti, Rec_true,Cyc_esti, Cyc_true,InTwo_esti, ... 
    InTwo_true, OutTwo_esti,OutTwo_true,Densi_esti,Densi_true] = evlua(X_output, X_observe, X_observe_full,Omiga)
%UNTITLED 
% MSE is evaluated at genuinely missing entries against the noiseless truth.
% R2 describes goodness of fit over the complete noiseless tensor.
missing = 1-double(Omiga);
error_full = double(X_output-X_observe_full);
n_missing = sum(missing(:));
if n_missing > 0
    MSE = sum((error_full(:).*missing(:)).^2)/n_missing;
else
    MSE = NaN;
end

truth = double(X_observe_full);
truth_mean = mean(truth(:));
sst = sum((truth(:)-truth_mean).^2);
if sst > 0
    R2 = 1-sum(error_full(:).^2)/sst;
else
    R2 = NaN;
end
%% Reciprocity
alpha = 0.8;

X = double(X_output);
[n,m,T] = size(X);


inner = zeros(T,1);
for t = 1:T
    A = X(:,:,t);
    B = A .* A.';                    % B_ij = A_ij * A_ji
    inner(t) = sum(B(:)) - sum(diag(B));  
end
Rec_esti = sum( abs(inner) .^ alpha )/T;

X = double(X_observe_full);
[n,m,T] = size(X);

alpha = 0.8;
inner = zeros(T,1);
for t = 1:T
    A = X(:,:,t);
    B = A .* A.';                    % B_ij = A_ij * A_ji
    inner(t) = sum(B(:)) - sum(diag(B));  % 去掉 i=j 的对角线
end
Rec_true = sum( abs(inner) .^ alpha )/T;

%% ﻿Cyclic triads
alpha = 0.8;

X = double(X_output);
[n,m,T] = size(X);


for t = 1:T
    A  = X(:,:,t);
    A2 = A*A;
    A3 = A2*A;

    trA3   = trace(A3);
    diagA  = diag(A);
    diagA2 = diag(A2);

    at = trA3 - 3*sum(diagA2 .* diagA) + 2*sum(diagA.^3);
    inner(t) = at;
end

Cyc_esti = sum( abs(inner).^ alpha )/T;

X = double(X_observe_full);
[n,m,T] = size(X);

for t = 1:T
    A  = X(:,:,t);
    A2 = A*A;
    A3 = A2*A;

    trA3   = trace(A3);
    diagA  = diag(A);
    diagA2 = diag(A2);

    at = trA3 - 3*sum(diagA2 .* diagA) + 2*sum(diagA.^3);
    inner(t) = at;
end

Cyc_true = sum( abs(inner).^ alpha )/T;

%% In two-stars 
alpha = 0.8;

X = double(X_output);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);

    colsum  = sum(A,1)';           % 每列求和 -> n×1
    d       = diag(A);             % 对角
    s       = colsum - d;          % s_i = sum_{r!=i} A_{ri}

    colsum2 = sum(A.^2,1)';        % 每列平方和
    q       = colsum2 - d.^2;      % q_i = sum_{r!=i} A_{ri}^2

    at = sum( s.^2 - q );          % a_t
    inner(t) = at;
end

InTwo_esti = sum( abs(inner) .^ alpha )/T;

X = double(X_observe_full);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);

    colsum  = sum(A,1)';           % 每列求和 -> n×1
    d       = diag(A);             % 对角
    s       = colsum - d;          % s_i = sum_{r!=i} A_{ri}

    colsum2 = sum(A.^2,1)';        % 每列平方和
    q       = colsum2 - d.^2;      % q_i = sum_{r!=i} A_{ri}^2

    at = sum( s.^2 - q );          % a_t
    inner(t) = at;
end

InTwo_true = sum( abs(inner) .^ alpha )/T;
%% Our two-strs
alpha = 0.8;

X = double(X_output);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);

    rowsum  = sum(A,2);        % 每行求和 -> n×1
    d       = diag(A);         % 对角
    s       = rowsum - d;      % s_i = sum_{c!=i} A_{ic}

    rowsum2 = sum(A.^2,2);     % 每行平方和
    q       = rowsum2 - d.^2;  % q_i = sum_{c!=i} A_{ic}^2

    at = sum( s.^2 - q );      % a_t
    inner(t) = at;
end

OutTwo_esti = sum( abs(inner) .^ alpha )/T;

X = double(X_observe_full);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);

    rowsum  = sum(A,2);        
    d       = diag(A);         % 
    s       = rowsum - d;      % s_i = sum_{c!=i} A_{ic}

    rowsum2 = sum(A.^2,2);     % 
    q       = rowsum2 - d.^2;  % q_i = sum_{c!=i} A_{ic}^2

    at = sum( s.^2 - q );      % a_t
    inner(t) = at;
end

OutTwo_true = sum( abs(inner) .^ alpha )/T;

%% Density
alpha = 0.8;

X = double(X_output);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);
    inner(t) = sum(A(:)) - trace(A);   % sum_{i≠j} A_ij(t)
end

Densi_esti = sum( abs(inner) .^ alpha )/T;

X = double(X_observe_full);
[n,m,T] = size(X);

for t = 1:T
    A = X(:,:,t);
    inner(t) = sum(A(:)) - trace(A);   % sum_{i≠j} A_ij(t)
end

Densi_true = sum( abs(inner) .^ alpha )/T;




end
