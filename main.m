%% U1

x0 = 0.75;

err_central = zeros(1,8);
err_left    = zeros(1,8);
err_right   = zeros(1,8);
hs          = zeros(1,8);

true_du = exp(x0)*cos(exp(x0));
% up1 = u0 plus 1 timestep

% Approx derivata med central metod
for k = 1:8
    h = 2^(-k);
    hs(k) = h;
    
    u = @(x) sin(exp(x));
    u0  = u(x0);
    up1 = u(x0 + h);
    up2 = u(x0 + 2*h);
    un1 = u(x0 - h);
    un2 = u(x0 - 2*h);

    % central
    d_central = (up1 - un1)/(2*h);

    % frammåt
    d_left = (-up2 + 4*up1 - 3*u0)/(2*h);

    % bakåt
    d_right = (3*u0 - 4*un1 + un2)/(2*h);

    % fel
    err_central(k) = abs(d_central - true_du);
    err_left(k)    = abs(d_left    - true_du);
    err_right(k)   = abs(d_right   - true_du);
end

figure(1); 
hold on;

loglog(1./hs, err_central)
loglog(1./hs, err_left)
loglog(1./hs, err_right)

% referenskurva O(h^2)
loglog(1./hs, hs.^2)

xlabel('1/h')
ylabel('Error')
legend('Central (2)', 'Forward (6)', 'Backward (7)', '1/h')
grid on
title('Error vs O(h^2')


%% U2

N = 400;
h = (2*pi)/N;
d = 5;
f = @(x) d*cos(2*x);
u = @(x) (delta-(d/12))*cos(4*x) + (gamma/4)*sin(4*x) + (d/12)*cos(2*x);
H = 1/(h^2);
x = linspace(0, 2*pi, N)';

A = zeros(N, N);
b = zeros(N, 1);
q = zeros(N, 1);

gamma = 1;
delta = 1;



% adding finite element method pattern to A
% adding function values to q
for i = 1:N-1
    A(i+1, i) = H;
    A(i+1, i+1) = 16-2*H;
    if i~=(N-1)
        A(i+1, i+2) = H;
    end

    q(i+1) = f(h*i);
end

% adding randvilkor to A
A(1, 1) = -1;
A(1, 2) = 1;

% adding randvilkor to q
q(1) = q(1) + gamma*h;
q(N) = q(N) - delta*H;


% Lös ekvationssystemet A*w = q -> w = A\q
w = A\q;

figure(2)
hold on
plot(x, w)
plot(x, u(x))

figure(3)
plot(x, w-u(x))

%%CN function

function [x, t, U] = crankNicolson1D(D, f_fun, TL_fun, TR_fun, u0_fun, T, N, M)
% Crank–Nicolson för ut = D u_xx + f på [0,1] med Dirichlet-randvillkor

dx = 1 / N;
dt = T / M;

x = linspace(0, 1, N+1)';      % rumspunkter
t = linspace(0, T, M+1);       % tidspunkter

U = zeros(N+1, M+1);

% initialdata
for j = 1:(N+1)
    U(j,1) = u0_fun(x(j));
end
U(1,1)   = TL_fun(0);
U(end,1) = TR_fun(0);

lambda = D * dt / dx^2;

% bygg A och B för de inre punkterna (j = 1,...,N-1)
A = zeros(N-1, N-1);
B = zeros(N-1, N-1);

for j = 1:(N-1)
    A(j,j) = 1 + lambda;
    B(j,j) = 1 - lambda;
    if j > 1
        A(j,j-1) = -lambda/2;
        B(j,j-1) =  lambda/2;
    end
    if j < N-1
        A(j,j+1) = -lambda/2;
        B(j,j+1) =  lambda/2;
    end
end

% inre punkter vid t = 0
w = U(2:N, 1);

% tidsstegning
for n = 1:M
    t_n   = t(n);
    t_np1 = t(n+1);

    uL_n   = TL_fun(t_n);
    uR_n   = TR_fun(t_n);
    uL_np1 = TL_fun(t_np1);
    uR_np1 = TR_fun(t_np1);

    % f i inre punkter
    f_n   = zeros(N-1, 1);
    f_np1 = zeros(N-1, 1);
    for j = 1:(N-1)
        f_n(j)   = f_fun(x(j+1), t_n);
        f_np1(j) = f_fun(x(j+1), t_np1);
    end

    % högerled
    rhs = B*w + 0.5*dt*(f_n + f_np1);

    % Dirichlet-randvillkor (bidrag från x=0 och x=1)
    rhs(1)   = rhs(1)   + 0.5*lambda*(uL_n + uL_np1);
    rhs(end) = rhs(end) + 0.5*lambda*(uR_n + uR_np1);

    % lös linjärt system
    w = A \ rhs;

    % spara hela profilen vid tiden t_{n+1}
    U(1,   n+1) = uL_np1;
    U(end, n+1) = uR_np1;
    U(2:N, n+1) = w;
end
end


%% U3

% parametrar
D  = 1;          % värmeledningskoefficient
T  = 0.01;       % slut-tid
N  = 100;        % antal delintervall i x
M  = 1000;       % antal tidssteg (dt = T/M = 1e-5)

dx = 1/N;
dt = T/M;

x  = linspace(0, 1, N+1)';      % x_0 = 0,  x_N = 1
t  = linspace(0, T, M+1);       % t_0 = 0,  t_M = T

% funktioner f, randvillkor och begynnelsevillkor
f_fun  = @(x, t) 0;             % f(x,t) = 0
TL_fun = @(t) 0;                % u(0,t) = 0
TR_fun = @(t) 0;                % u(1,t) = 0
u0_fun = @(x) sin(5*pi*x);      % u(x,0) = sin(5*pi*x)

% initialdata vid t = 0
u = zeros(N+1, 1);
for j = 1:(N+1)
    u(j) = u0_fun(x(j));
end
u(1)   = TL_fun(0);
u(end) = TR_fun(0);

% enkel stabilitetskontroll
r = D*dt/(dx^2);
if r > 0.5
    warning('D*dt/dx^2 = %f > 1/2, schemat kan vara instabilt.', r)
end

% tidsstepping med Euler framåt (explicit)
for n = 1:M
    t_n   = t(n);
    t_np1 = t(n+1);

    % randvillkor vid tiden t_n
    u(1)   = TL_fun(t_n);
    u(end) = TR_fun(t_n);

    u_new = u;

    % inre punkter j = 2,...,N
    for j = 2:N
        lap = u(j+1) - 2*u(j) + u(j-1);    % andraderivata med central differens
        u_new(j) = u(j) + dt*( D*lap/(dx^2) + f_fun(x(j), t_n) );
    end

    % randvillkor vid tiden t_{n+1}
    u_new(1)   = TL_fun(t_np1);
    u_new(end) = TR_fun(t_np1);

    u = u_new;
end

% analytisk lösning
% för f=0, D=1, u0(x)=sin(5*pi*x)
u_exact = exp(-(5*pi)^2 * D * T) * sin(5*pi*x);

% fel vid tiden T
err = u - u_exact;

figure
plot(x, err)
xlabel('x')
ylabel('u_{num} - u_{exact}')
title('U3(c): fel vid T = 0.01')
grid on

%% U4c

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
