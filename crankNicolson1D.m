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
