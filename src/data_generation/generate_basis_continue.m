function output = generate_basis_continue(K)
%GENERATE_BASIS_CONTINUE Used to generate a set of Fourier orthogonal basis functions for observations on [-1, 1]
%input:K The number of orthogonal basis functions to generate
%output:A 1×K cell array, each element is a function handle;
%       the k-th element returns cos(pi*(k-1)*x)

output = cell(1,K);
output{1} = @(x)sqrt(2)/2;
fourier = @(x,k)cos(k*pi*x);
for k = 2:K
   output{k} = @(x)fourier(x,k-1); 
end
end

