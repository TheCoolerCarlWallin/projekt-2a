%% U4(f)


uanafcn  = @(x,t) sin(exp(t).*x);
uinitfcn = @(x)   uanafcn(x,0);
ffcn     = @(x,t) exp(t).*x.*cos(exp(t).*x) + exp(2*t).*sin(exp(t).*x);
TLfcn    = @(t)   0;
TRfcn    = @(t)   sin(exp(t));

D     = 1;
N     = 200;         
dx    = 1/N;
dt0   = 0.5 * dx;     
Tvals = [1 2 4 8];    

for m = 1:length(Tvals)
    T  = Tvals(m);
    M  = round(T / dt0);
    dt = T / M;      

    [x, t, U] = crankNicolson1D(D, ffcn, TLfcn, TRfcn, uinitfcn, T, N, M);
    u_exact   = uanafcn(x, T);

    figure
    plot(x, u_exact, 'k--', x, U(:,end), 'b-');
    legend('Exakt','Numerisk','Location','Best');
    xlabel('x');
    ylabel('u(x,t)');
    title(sprintf('U4(f): t = %.1f  (dt = %.3e)', T, dt));
    grid on
end
