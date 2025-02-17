clear all

k = 10;

N = 2 * k + 1;
N_num = 9999;

sample = zeros(1, N-1);
for j = 1:N-1;

    sample(j) = H_fun((j - 1) * (1/N));

end


sample(2) = -1i / (2 * k + 1);
sample(k) = -k * 1i / (2 * k + 1);
sample(k + 1) = (k + 1) * 1i / (2 * k + 1);
sample(2 * k) = 2i * k / (2 * k + 1);

r1n = ifft(sample);
rn = [r1n(ceil(N / 2) : end) , r1n(1 : floor(N / 2) + 1)];


F = linspace(0.0 , 1.0 , N_num);

RF = zeros(1 , N_num);

for j = 1 : N_num;
    a = 0;
    for n = -k:k
        a = a + rn(n + k + 1) * exp(-1i * 2 * pi * F(j) * n);
    end
    RF(j) = imag(a);
end

figure;
plot(F, RF, 'LineWidth', 2);
hold on;
plot(F, imag(arrayfun(@H_fun, F)), 'LineWidth', 2);
hold off;
title('Frequency Response');
legend('R(F)', 'H_d(F)');
xlabel('Frequency (F)');
ylabel('Magnitude');
xlim([-0.1 1.1])
ylim([-1.5 1.5])
grid on;

figure;
stem_x = [-k : k];
stem( stem_x, rn, 'LineWidth', 2);
txt = ['N = ' , num2str(length(rn))];
text(-2.5 , 0.85 , txt,'FontSize' , 14)
title('Impulse Response r[n]');
legend('r[n]');
xlabel('Sample index (n)');
ylabel('Magnitude');
xlim([-k-1 k+1])
ylim([-1 1])
grid on;

function Hd = H_fun(x)
    if x == 0
        Hd = 0;
    elseif x > 0 && x <= 0.5
        Hd = -1i;
    elseif x > 0.5 && x <= 1
        Hd = 1i;
    else
        Hd = NaN; % Handle out-of-range values
    end
end
