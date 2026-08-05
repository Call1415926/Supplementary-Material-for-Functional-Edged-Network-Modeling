function [output, U_true,B] = generate_data_true(n,m,s1,s2,K,N,t_start, t_end, seq)
%GENERATE_DATA Used to generate a set of initial functions to be observed
%input:n, m The first and second dimensions of X; the third functional dimension uses Fourier basis and is fixed as [-1, 1]
%     :s1 s2 K The core dimensions of the factor matrices in the decomposition
%     :N The maximum allowed value when generating the Tucker core, taking values in [0, N]
%output:An n×m cell array; each element is a function handle
%      :Psi The basis for the first mode of the data
%      :Phi The basis for the second mode of the data
%      :Theta_dis The full discrete observations of the basis for the third mode of the data
%      :B The core of the generated data, a 3D tensor

Psi = generate_basis_discrete(s1,n);
Phi = Psi;%generate_basis_discrete(s2,m); 
Theta = generate_basis_continue(K);
B = generate_core(s1,s2,K,N);
X_discrete = ones(n,m,K);
for k = 1:K
   X_discrete(:,:,k) = Psi * B(:,:,k) * Phi'; 
end
seqtime = t_start:seq:t_end;
myseq = length(seqtime);
Theta_dis = ones(myseq, K);
for k = 1:K
   for j = 1:myseq
      Theta_dis(j, k) = Theta{k}(seqtime(j)) ;
   end
end

U_true = {Psi,Phi, Theta_dis};

X_continue = cell(n,m);
for i = 1:n
    for j = 1:m
        vector_K = reshape(X_discrete(i,j,:),[],1);
        X_continue{i,j} = @(x)generate_data_true_funhandle(x,vector_K,Theta);
    end
end
B = tensor(B);
output = X_continue;
end

