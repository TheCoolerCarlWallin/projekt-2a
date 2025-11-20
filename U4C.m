%% U4

% Exakta funktioner enligt (18)–(19)
uanafcn  = @(x,t) sin(exp(t) .* x);
uinitfcn = @(x)   uanafcn(x, 0);
ffcn     = @(x,t) exp(t) .* x .* cos(exp(t) .* x) + exp(2*t) .* sin(exp(t) .* x);
TLfcn    = @(t)   0;
TRfcn    = @(t)   sin(exp(t));

D = 1;
T = 0.1;
N = 100;
M = 100;          

[x, t, U] = crankNicolson1D(D, ffcn, TLfcn, TRfcn, uinitfcn, T, N, M);

% analytisk lösning vid T
u_exact = uanafcn(x, T);

err = U(:, end) - u_exact;

fprintf('Maxfel vid T = %.3f: %.3e\n', T, max(abs(err)));

figure
plot(x, err)
xlabel('x')
ylabel('U_{num} - u_{exact}')
title('U4(c): fel vid T = 0.1')
grid on
