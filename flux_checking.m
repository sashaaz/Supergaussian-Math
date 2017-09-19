%n = 1;
d_fwhm = 100e-6;
%fun = {};

%for n = 1:100
%    fun{n} = @(x) x/(d_fwhm)^2 .* exp(-(2*x/d_fwhm).^(2*n) * log(2));
%    fun{n} = integral(fun{n}, 0, Inf);
%end
%fun_values = cell2mat(fun);

%for n = 1:100
%    A0(n) = sqrt(gamma(1/n)/(n*log(2)^(1/n)));
%end

%plot(A0)
%hold on
%plot(fun_values + 1 - .13);
%legend('Integral of x * exp(-(2*x/d_{fwhm})^{2n} * log(2))', 'Dependence of 1/A0');
%xlabel('n');


x = d_fwhm*(0.01 : 0.01 : 100);
dx = x(2) - x(1)
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

