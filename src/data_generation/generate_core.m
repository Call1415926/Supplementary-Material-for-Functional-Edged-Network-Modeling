function output = generate_core(s1,s2,K,N)
%GENERATE_CORE Generate a 3D tensor whose slices are as full-rank as possible, 
%              with each element taking integer values from 0 to N; for fixed K, 
%              each slice along the third mode is full-rank
%input:s1 s2 K The three dimensions of the tensor
%output:The generated tensor
output = ones(s1,s2,K);
for k = 1:K
   p = 0;
   while p == 0
       randmatrix = ceil(rand(s1,s2)*(N-2) + 1);
       if rank(randmatrix) == min(s1,s2)
          p = 1;
          output(:,:,k) = randmatrix;
       end
   end
end
end

