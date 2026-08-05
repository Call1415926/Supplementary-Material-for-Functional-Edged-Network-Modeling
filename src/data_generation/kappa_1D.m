function output = kappa_1D(x)
output = 0;
if abs(x) <= 1
   output = 0.75 * (1 - x^2); 
end
end

