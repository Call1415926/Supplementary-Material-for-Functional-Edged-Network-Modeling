function [X_observe, Omiga, seqtime,X_observe_full, X_observe_full_nornd] = generate_data_observe(X_true,N_min,N_max,t_start,t_end,seq,sigma)
%GENERATE_DATA_OBSERVE Function to generate observations according to the true data distribution
%input:X_true An n×m cell array, each element is a function handle representing the corresponding function
%     :N_min, N_max The minimum and maximum number of observation points in each dimension;
%                   note that the maximum number of observations must be less than
%                   the number of subintervals after partitioning
%     :t_start, t_end The function interval for observations; here [-1, 1] is used
%                     as the start and end interval for Fourier basis observations
%     :seq The length of each small subinterval, i.e., the spacing between grid points
%     :sigma The standard deviation of the observation noise
%output:X_observe The observed data, a 3D tensor where unobserved entries are filled with 0;
%                 the first two modes match the size of X_true, and the size of the third mode
%                 is determined jointly by t_start, t_end, and seq
%       X_observe_full The fully observed data with random noise, a 3D tensor
%       X_observe_full_nornd The fully observed data without random noise, a 3D tensor
%       Omiga A 0/1 tensor with the same size as X_observe; observed entries are 1, others are 0
%       seqtime The true observation time points, a vector whose length is determined jointly by
%               t_start, t_end, and seq

myseq = t_start:seq:t_end;
mylength = length(myseq);
[n,m] = size(X_true);
X_observe_full = ones(n, m, mylength);
Omiga = zeros(n, m, mylength);
for i = 1:n
   for j = 1:m
       N_observe = unidrnd(N_max - N_min + 1) + N_min - 1; %生成随机观测的数量
       temp = randperm(mylength);
       Omiga(i, j, temp(1:N_observe)) = 1;
      for k = 1:mylength
         X_observe_full(i, j, k) = X_true{i, j}(myseq(k)) + normrnd(0, sigma);
         X_observe_full_nornd(i, j, k) = X_true{i, j}(myseq(k));
      end
   end
end
Omiga = tensor(Omiga); %部分观测映射张量
X_observe_full = tensor(X_observe_full);
X_observe_full_nornd = tensor(X_observe_full_nornd);
X_observe = times(X_observe_full, Omiga); 
seqtime = myseq;
end

