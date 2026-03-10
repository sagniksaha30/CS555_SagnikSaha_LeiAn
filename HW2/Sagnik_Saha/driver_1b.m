Ts = 0.1:0.01:2;      % T values
%figure;
figure()
%s(t) array
st = zeros(length(Ts),1);
for j = 1:length(Ts)
    T = Ts(j);
    run('bd.m')
    %Take derivative
    w = D*u;
    st(j) = max(abs(w(:)));
end
plot(Ts,st, 'LineWidth', 1);
title('$\nu = 0.003183, dt=0.0005$','Interpreter', 'latex','FontSize',8);
axis square;
grid on;
xlabel('t');
ylabel('s(t)');
%ylim([0.0 1.0]);

[stmax, ind] = max(st);
disp(['s_* = ', num2str(stmax)]);
disp(['t_* = ', num2str(Ts(ind))]);

%legend('Location', 'eastoutside'); 
%exportgraphics(gcf,'/Users/sagniksaha/Desktop/CS555/hw/tex/hw2/figs/1a.pdf','ContentType','vector');
%print('/Users/sagniksaha/Desktop/CS555/hw/tex/hw2/figs/1b.pdf','-dpdf','-r300')
pause;