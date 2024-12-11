%% U1: Spektral differetiering
clc; clear;

%% Parametrar och funktioner

alpha = 400;
f = @(x) exp(-alpha * (x - pi/5).^2);                           % f(x)
df = @(x) -2 * alpha * (x - pi/5) .* f(x);                      % f'(x)
d2f = @(x) 2 * alpha * f(x) .* (2 * alpha * (x - pi/5).^2 - 1); % f''(x)
N_values = 2.^(4:8);                                            % N = 2^4, 2^5, ..., 2^8

%% Main code

errors_df = zeros(size(N_values));                              % Fel för första derivatan
errors_d2f = zeros(size(N_values));                             % Fel för andra derivatan

for i = 1:length(N_values)
    
    N = N_values(i);                                            % Antal punkter
    x = (0:N-1)/N;                                              % Diskretisering av intervallet [0, 1]
    dx = 1/N;                                                   % Steglängd

    % Värden av funktionen f vid de diskreta punkterna
    f_values = f(x)';

    % Spektralmetod för att approximera första derivatan
    f_hat = fft(f_values);                                      % Fourierkoefficienter
    k = [0:N/2 -N/2+1:-1]';                                     % Frekvensvektor för jämnt N
    df_hat = 2i * pi * k .* f_hat;                              % Fourierkoefficienter för f'(x)
    df_approx = real(ifft(df_hat));                             % Tillbaka till fysiskt rum

    % Spektralmetod för att approximera andra derivatan
    d2f_hat = -(2 * pi * k).^2 .* f_hat;                        % Fourierkoefficienter för f''(x)
    d2f_approx = real(ifft(d2f_hat));                           % Tillbaka till fysiskt rum

    % Exakta värden av derivator
    df_exact = df(x)';
    d2f_exact = d2f(x)';

    % Beräkna fel i diskret 2-norm
    errors_df(i) = sqrt(sum((df_approx - df_exact).^2) / N);
    errors_d2f(i) = sqrt(sum((d2f_approx - d2f_exact).^2) / N);
end

% Plot av fel i log-log
figure;
loglog(N_values, errors_df, 'r--o'); hold on;
loglog(N_values, errors_d2f, 'b--s');
legend({'$\|Df - f''\|_2$', '$\|D^2f - f''''\|_2$'}, 'Interpreter', 'latex');
xlabel('$N$', 'Interpreter', 'latex');
ylabel('Fel i diskret 2-norm', 'Interpreter', 'latex');
title('Fel i approximation av derivator (log-log)');
grid on;

