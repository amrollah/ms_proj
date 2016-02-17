function y = vml_cc2pow(x)
y = [1-x 1+(1-exp(x))/(exp(1)-1) (1-tanh(4*(x-0.5))/tanh(2))/2]';
y = y(:);

