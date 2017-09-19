n = 1;
d_fwhm = 100e-6;
fun = {};

for n = 1:100
    fun{n} = @(x) x/(d_fwhm)^2 .* exp(-(2*x/d_fwhm).^(2*n) * log(2));
    fun{n} = integral(fun{n}, 0, Inf);
end
fun_values = cell2mat(fun);

for n = 1:100
    A0(n) = sqrt(gamma(1/n)/(n*log(2)^(1/n)));
end

plot(A0)
hold on
plot(fun_values + 1 - .13);
legend('Integral of x * exp(-(2*x/d_{fwhm})^{2n} * log(2))', 'Dependence of 1/A0');
xlabel('n');

