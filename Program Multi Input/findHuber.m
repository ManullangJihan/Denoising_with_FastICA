%% Fungsi Huber
function [G, g, g_prime] = findHuber(x) 

% Fungsi Estimasi Ng_Entropy
a = 1;

% Fungsi Huber
G = log(cosh(x));

% Turunan Pertama Fungsi Huber
g = tanh(a*x);

% Turunan Kedua Fungsi Huber
g_prime = a*(1-(tanh(a*x)).^2);

end