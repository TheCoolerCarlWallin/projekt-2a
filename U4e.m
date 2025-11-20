%% U4(e) – konvergensstudie

uanafcn  = @(x,t) sin(exp(t).*x);
uinitfcn = @(x)   uanafcn(x,0);
ffcn     = @(x,t) exp(t).*x.*cos(exp(t).*x) + exp(2*t).*sin(exp(t).*x);
TLfcn    = @(t)   0;
TRfcn    = @(t)   sin(exp(t));

D  = 1;
T  = 1;                        
Ns = [25 50 100 200 400];      

errs  = zeros(size(Ns));
dxs   = zeros(size(Ns));
dts   = zeros(size(Ns));

c_dt = 0.2;                   

for k = 1:length(Ns)
    N  = Ns(k);
    dx = 1/N;
    dxs(k) = dx;
    
    dt_guess = c_dt * dx;
    M  = round(T / dt_guess);  
    dt = T / M;
    dts(k) = dt;

    [x, t, U] = crankNicolson1D(D, ffcn, TLfcn, TRfcn, uinitfcn, T, N, M);

    u_exact = uanafcn(x, T);

    % diskret 2-norm över inre punkter (2..N)
    r       = U(2:N,end) - u_exact(2:N);
    errs(k) = sqrt( (1/(N-1)) * sum(abs(r).^2) );
end

% uppmätt konvergensordning
orders = NaN(size(Ns));
for k = 2:length(Ns)
    orders(k) = log(errs(k-1)/errs(k)) / log(dxs(k-1)/dxs(k));
end

fprintf('\nU4(e): konvergensstudie vid T = %.1f\n', T);
fprintf('   N        dx          dt          fel(2-norm)    ordning\n');
for k = 1:length(Ns)
    if k == 1
        ord_str = '   -';
    else
        ord_str = sprintf('%10.3f', orders(k));
    end
    fprintf('%4d   %9.4e  %9.4e  %13.4e  %s\n', Ns(k), dxs(k), dts(k), errs(k), ord_str);
end
