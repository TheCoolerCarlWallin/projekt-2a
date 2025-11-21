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

