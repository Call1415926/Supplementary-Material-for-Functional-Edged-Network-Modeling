function output = dif_mat(n)
%DIF_MAT Generate the difference matrix
%input:n The number of rows of the square difference matrix
%output:The resulting n×n difference matrix
D = zeros(n-1, n);
for i = 1:(n - 1)
   D(i, i) = 1;
   D(i, i + 1) = -1;
end
output = D'*D;
end

