d_fwhm = 100e-6;

%% Analytical method of checking.
fun = {};

for n = 1:100
   fun{n} = @(x) x/(d_fwhm)^2 .* exp(-(2*x/d_fwhm).^(2*n) * log(2));
   fun{n} = integral(fun{n}, 0, Inf);
end
fun_values = cell2mat(fun);

for n = 1:100
   A0(n) = sqrt(n*log(2)^(1/n)/(gamma(1/n)));
end
P0 = abs(A0).^2 .* fun_values;
plot(n, P0);
xlabel('n');
title('P0 - i.e. |A0|^2 * integral (should be constant)');

figure
plot(0.125*abs(A0).^(-2))
hold on
plot(fun_values)
title('Curve of n-dependence for .125/|A0|^2 and the integral - should be the same');
legend('.125/|A0|^2', 'integral');


%% Numerical Method of checking
clear all

x = d_fwhm*(0.01 : 0.01 : 100);
dx = x(2) - x(1);
for n = 1:100
    r_0 = d_fwhm / 2 / (log(2)).^(1/2/n);
    f_x = exp(-(x / r_0).^(2*n)) .* x;
    int_f(n) = trapz(f_x) * dx;

    A_02(n) = 4 * log(2)^(1/n) * n / d_fwhm^2 / gamma(1/n);
end

n = 1:100;
figure;
plot(n, int_f); title('numerically integrated beam power when d_FWHM is constant'); xlabel('n')
figure;
plot(n, A_02);  title('analytically computed |A_0|^2 to keep P_0 and d_FWHM constant'); xlabel('n')
figure;
plot(n, int_f .* A_02); title('P_0 = |A_0|^2 * integral (should be constant)'); xlabel('n')


r_0 = d_fwhm(1) / 2 / (log(2)).^(1/2)
x = r_0*(0.01 : 0.01 : 100);
dx = x(2) - x(1)
for n = 1:100
    f_x = exp(-(x / r_0).^(2*n)) .* x;
    int_f(n) = trapz(f_x) * dx;

    A_02(n) = n / r_0^2 / gamma(1/n);
end

n = 1:100;
figure;
plot(n, int_f); title('numerically integrated beam power when r_0 is constant'); xlabel('n')
figure;
plot(n, A_02);  title('analytically computed |A_0|^2 to keep P_0 and r_0 constant'); xlabel('n')
figure;
plot(n, int_f .* A_02); title('P_0 = |A_0|^2 * integral (should be constant)'); xlabel('n')

