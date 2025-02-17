clear all
N = 1200;
M = 2;
L = [1:500];
y = N ./ L * 3 .* (L+M-1) .* (log2(L + M - 1) + 1);
[C0 , L0] = min(y);