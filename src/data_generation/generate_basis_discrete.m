function output = generate_basis_discrete(s,n)
%Generate a set of discrete orthonormal basis vectors
%input:s The number of basis vectors to generate
%      n The dimension of each basis vector
%output:An n×s real matrix, where each column is one basis vector


A = rand(n,n);
B = orth(A);
output = B(:,1:s);
end