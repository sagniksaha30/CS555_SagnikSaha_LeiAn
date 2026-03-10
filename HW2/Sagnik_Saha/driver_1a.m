Ts = 0.1:0.1:2;      % T values
%figure;
figure('Position',[100, 100, 800, 1000])
hold on;
for k = 1:length(Ts)
    T = Ts(k);
    run('bd.m')
    
    % random color
    c = rand(1,3);  % RGB triplet with values between 0 and 1
    
    % plot with random color
    plot(xb, ufinal, 'Color', c, 'LineWidth', 1,'DisplayName', sprintf('T = %.1f', T));
end

title('Burgers Equation, IMEX BDF3/EXT3, $\nu = 0.01, dt=0.0005$','Interpreter', 'latex','FontSize',8);
axis square;
grid on;
xlabel('x');
ylabel('u(x,T)');
ylim([0.0 1.0]);

legend('Location', 'eastoutside'); 
%exportgraphics(gcf,'/Users/sagniksaha/Desktop/CS555/hw/tex/hw2/figs/1a.pdf','ContentType','vector');
print('/Users/sagniksaha/Desktop/CS555/hw/tex/hw2/figs/1a.pdf','-dpdf','-r300')
pause;