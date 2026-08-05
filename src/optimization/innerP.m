function output = innerP(A,B)
%INNERP Used to compute the inner product of two tensors
%input:A, B Two tensors of the same size
%output:A scalar equal to the inner product of A and B

C = times(A, B);
C = double(C);
output = sum(C(:));


end

