function output = generate_data_true_funhandle(x,vector_K, cell_fun)
%GENERATE_DATA_TRUE_FUNHANDLE A helper function used inside generate_data_true, returning a function handle
%input:x The time parameter to be evaluated
%     :vector_K A data vector of length K
%     :cell_fun A cell array of length K, each element being a function handle

output = 0;
K = length(vector_K);
for k = 1:K
   output = output + vector_K(k) * cell_fun{k}(x);
end
output = output + mu_true(x);
end

