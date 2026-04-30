Ex_values = [20, 20*2, 20*4, 20*8, 20*16,20*32,20*64];
results= zeros(size(Ex_values));

for j = 1:length(Ex_values)
    Ex = Ex_values(j);
    run('taper_flat_a.m');
    results(j) = mpe;
end

loglog(Ex_values, results, 'o-');
hold on;
loglog(Ex_values, (10^5)*Ex_values.^(-2), 'r--', 'LineWidth', 1);
xlabel('Ex');
legend('Max pointwise error', '$y = 10^{5}/E_x^2$', 'Interpreter', 'latex');
grid on
saveas(gcf, 'p1b_convergence.png');
pause()
