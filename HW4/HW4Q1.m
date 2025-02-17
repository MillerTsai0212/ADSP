clear all
close all

A = imread('parrots.jpg');
A = rgb2gray(A);


B = imread('parrots.jpg');
B = rgb2gray(B);


c1 = 1 / (255)^(1/2);
c2 = 1 / (255)^(1/2);


[B, ssim] = SSIM(A, B, c1, c2);




L=255;
function [B, ssim] = SSIM(A, B, c1, c2)
    
    L = 255;
    
   
    B = B * 0.5 + 255.5 * 0.5;
    
    
    B = uint8(B);
    
   
    [M, N] = size(A);
    
    
    mu_x = mean(A(:));
    mu_y = mean(B(:));
    
    
    sig_x = sum((double(A(:)) - mu_x).^2) * (1 / (M * N));
    sig_y = sum((double(B(:)) - mu_y).^2) * (1 / (M * N));
    
    
    sig_xy = sum((double(A(:)) - mu_x) .* (double(B(:)) - mu_y)) * (1 / (M * N));
    
   
    ssim = (2 * mu_x * mu_y + (c1 * L)^2) * ...
                 (2 * sig_xy + (c2 * L)^2) / ...
                 ((mu_x^2 + mu_y^2 + (c1 * L)^2) * ...
                 (sig_x + sig_y + (c2 * L)^2));
    
    
    fprintf('SSIM = %.4f\n', ssim);
    
   
    figure;
    subplot(1, 2, 1), imshow(A), title('Image A');
    subplot(1, 2, 2), imshow(B), title('Image B');
end

