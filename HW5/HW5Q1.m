clear all
x = [1, 2, 3, 4]; 
y = [1, 2, 3, 4];

[Fx, Fy] = fftreal(x, y);

disp('訊號號 x 的 FFT:');
disp(Fx);

disp('訊號 y 的 FFT:');
disp(Fy);
function [Fx, Fy] = fftreal(x, y)
    % 確保 x 和 y 是列向量
    x = x(:);
    y = y(:);

    % 將 x 和 y 合併成一個複數信號 z
    z = x + 1i * y;

    % 計算複數信號 z 的 FFT
    Z = fft(z);

    % 提取 x 和 y 的 FFT
    N = length(x);
    Z_conj_flip = conj(flipud(Z));
    
    % 分別計算 x 和 y 的 FFT
    Fx = 0.5 * (Z + Z_conj_flip);
    Fy = -1i * 0.5 * (Z - Z_conj_flip);

    % 取實部，確保結果為實數
    Fx = real(Fx);
    Fy = real(Fy);
end