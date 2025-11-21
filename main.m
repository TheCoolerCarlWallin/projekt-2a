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
