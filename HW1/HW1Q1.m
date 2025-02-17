clear all

N = 17;
delta = 1e-4;
t1 = 0: delta: 0.5;
t2 = 0: N-1;
HD_PASS = [0, 0.25];
F_m= [0, 0.05, 0.1, 0.15, 0.19, 0.26, 0.3, 0.35, 0.4, 0.45];
k = (N-1)/2;
N_EXP = k + 2;
m_err = 99;

while true
    A = [];
    for m = 0:k+1
        row = [];
        for n = 0:k+1
            if n == 0
                row = [row, 1];
            elseif n == k+1
                row = [row, ((-1)^m) / W(F_m(m+1))];
            else
                row = [row, cos(2*n*pi*F_m(m+1))];
            end
        end
        A = [A; row];
    end
    
    A_inv = inv(A);
    
    b = [];
    for m = 0:k+1
        b = [b, Hd(F_m(m+1))];
    end
    
    b = b';
    
    S = A_inv * b;
    
    % Step3 find error's local maximum
    ex_n = [];
    max_err = -1;
    F_ll = [];
    F_l = [];
    for i = 1:(0.5 / delta) + 1 % TODO add zero to the begin and the end
        if i == 1
            F_ll = 0;
            F_l = err(0 * delta, S);
        end
        
        if i == (0.5 / delta) + 1
            F = 0;
        else
            F = err(i * delta, S);
        end
        
        % check if F_l is a minimal or maximum
        if (F_l - F > 0 && F_l - F_ll > 0) || ...
           (F_l - F < 0 && F_l - F_ll < 0) % Maximum or Minimum
            ex_n = [ex_n, (i - 1) * delta];
            if max_err < abs(F_l)
                max_err = abs(F_l);
            end
        end
        
        % update F_l and F_ll
        F_ll = F_l;
        F_l = F;
    end
    
    fprintf('max_err = %f\n', max_err);
    
    if 0 <= m_err - max_err && m_err - max_err <= delta
        break;
    end
    
    % Update m_err and F[m]
    m_err = max_err;
    F_m = ex_n(1:k + 2);
end

h = calculate_h(N, k, S);



%% Plot
% Plot RF(t) using S values
subplot(211)
RF_values = zeros(1, length(t1));
for i = 1:length(t1)
    RF_values(i) = RF(t1(i), S);
end
plot(t1, RF_values, 'LineWidth', 2);
hold on;

% Plot Hd(t)
Hd_values = zeros(1, length(t1));
for i = 1:length(t1)
    Hd_values(i) = Hd(t1(i));
end
plot(t1, Hd_values, 'LineWidth', 2);
legend('FIR Filter', 'Designed Filter');
title('Frequency Response');
grid on;
hold off;


subplot(212)
stem(t2, h, 'LineWidth', 2);
grid on;
legend('h[n]');
title('Impluse Response');



%% Function
function output = W(x)
    WEIGHT = [0 0.2 10; 0.25 0.5 6];  % 定義權重矩陣
    pass_band_start = WEIGHT(1, 1);
    pass_band_end = WEIGHT(1, 2);
    pass_band_weight = WEIGHT(1, 3);
    
    stop_band_start = WEIGHT(2, 1);
    stop_band_end = WEIGHT(2, 2);
    stop_band_weight = WEIGHT(2, 3);

    if x >= pass_band_start && x <= pass_band_end
        output = pass_band_weight;
    elseif x >= stop_band_start && x <= stop_band_end
        output = stop_band_weight;
    else
        output = 0;
    end
end
function output = Hd(x, HD_PASS)
    HD_PASS = [0 0.225; 0.275 1];
    pass_band_start = HD_PASS(1, 1);
    pass_band_end = HD_PASS(1, 2);
    stop_band_start = HD_PASS(2, 1);
    stop_band_end = HD_PASS(2, 2);

    if x >= pass_band_start && x <= pass_band_end
        output = 1;
    else
        output = 0;
    end
end
function output = err(F, S, HD_PASS, WEIGHT)
    output = (RF(F, S) - Hd(F)) * W(F);
end
function output = RF(F, S)
    output = 0;
    k = 8;
    for n = 0:k
        ans_n = S(n+1) * cos(2 * pi * n * F);
        output = output + ans_n;
    end
end
function h = calculate_h(N, k, S)
    h = zeros(1, N);
    for i = 1:N
        if i < k+1
            h(i) = S(k-i+2) / 2;
        elseif i == k+1
            h(i) = S(1);
        else  % i > k+1
            h(i) = S(i-k) / 2;
        end
    end
end
