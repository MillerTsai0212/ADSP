clear all
close all

A =imread('parrots.jpg');
B = C420(A);
imwrite(uint8(B),'parrots_420.jpg','jpg');

function B=C420(A)
A = double(A);
[m,n,h] = size(A);

Yuv_coe = [0.299, 0.587, 0.114;
           -0.169, -0.331, 0.500;
           0.500, -0.419, -0.081];

Yuv = zeros(m,n,h);
sample = zeros(m,n,h);
B = zeros(m,n,h);

for i = 1 : m
    for k = 1 : n
        Yuv(i, k, :) = Yuv_coe * squeeze(A(i, k, :));
    end
end


for i = 1:m
    for j = 1:n
        sample(i, j, 1) = Yuv(i, j, 1);
        if mod(i, 2) == 1 && mod(j, 2) == 1
            sample(i, j, 2) = Yuv(i, j, 2);
        elseif mod(i, 2) == 0 && mod(j, 2) == 1
            sample(i, j, 3) = Yuv(i, j, 3);
        end
    end
end


for j = 1:2:n
    sample(1, j, 3) = Yuv(2, j, 3);
    sample(m, j, 2) = Yuv(m-1, j, 2);
end


for i = 2:2:size(sample, 1)-2
    for j = 1:2:size(sample, 2)
        temp = (sample(i-1, j, 2) + sample(i+1, j, 2)) / 2;
        sample(i, j, 2) = temp;
    end
end


for i = 3:2:size(sample, 1)
    for j = 1:2:size(sample, 2)
        temp = (sample(i-1, j, 3) + sample(i+1, j, 3)) / 2;
        sample(i, j, 3) = temp;
    end
end


for i = 2:2:size(sample, 2)-2
    for j = 1:size(sample, 1)
        for k = 2:3
            temp = (sample(j, i-1, k) + sample(j, i+1, k)) / 2;
            sample(j, i, k) = temp;
        end
    end
end


for n = 3:size(sample, 2)-1
    sample(:, n, 2:3) = (sample(:, n-1, 2:3) + sample(:, n+1, 2:3)) / 2;
end


for i = 1:m
    for k = 1:n
        
        sample_element = [sample(i, k, 1); sample(i, k, 2); sample(i, k, 3)];
       
        transformed_element = inv(Yuv_coe) * sample_element;
        
        B(i, k, :) = transformed_element';
    end
end

B = uint8(B);

end