clc; clear all; close all;

% ====== PARAMETER ======
c  = 1;
kd = 0.3;
l  = 1;
h  = 1;
dx = 0.05;
dt = 0.008;

% ====== GRID ======
x  = 0:dx:l;
t  = 0:dt:1.5;
Nx = length(x);
Nt = length(t);

% ====== KOEFISIEN A dan r ======
A  = 0.5*c^2 + 2*c^2/l;
r  = A*dt^2/dx^2;
cw = sqrt(A);
fprintf('A = %.4f\n', A);
fprintf('r = %.4f\n', r);

% ====== KOEFISIEN BTCS (dari materi) ======
Da = 2*l*dx^2 + 2*cw^2*l*dt^2 + 8*cw^2*dt^2 + 2*kd*l*dx^2*dt;
Db = l*dx^2   + cw^2*l*dt^2   + 4*cw^2*dt^2 + kd*l*dx^2*dt;

aa = (cw^2*l*dt^2 + 4*cw^2*dt^2) / Da;
bb = (2*l*dx^2    + kd*l*dx^2*dt) / Db;
cc =  l*dx^2 / Db;

fprintf('a = %.5f\n', aa);
fprintf('b = %.5f\n', bb);
fprintf('c = %.5f\n', cc);

% ====== KONDISI AWAL ======
u = zeros(Nx, Nt);
for j = 1:Nx
    if x(j) <= 0.5*l
        u(j,1) = 2*h*x(j)/l;
    else
        u(j,1) = 2*h - 2*h*x(j)/l;
    end
end

% ====== KONDISI BATAS ======
for n = 1:Nt
    u(1,n)  = 0;
    u(Nx,n) = 0;
end

% ====== INISIALISASI n=1 ======
for j = 2:Nx-1
    u(j,2) = u(j,1) + 0.5*r*(u(j+1,1) - 2*u(j,1) + u(j-1,1));
end

% ====== MATRIKS TRIDIAGONAL A ======
M     = Nx - 2;
A_tri = diag(ones(M,1)) ...
      - aa*diag(ones(M-1,1), 1) ...
      - aa*diag(ones(M-1,1),-1);

% ====== ITERASI IMPLISIT BTCS ======
for n = 2:Nt-1
    B = bb*u(2:Nx-1, n) - cc*u(2:Nx-1, n-1);
    u(2:Nx-1, n+1) = A_tri \ B;
end

% ====== PLOT 2D ======
figure(1);
hold on; grid on;
for n = 1:Nt
    plot(x, u(:,n), 'r-', 'LineWidth', 1);
end
xlabel('x (m)');
ylabel('u(x,t)');
title('Vibrasi Sasando — Implisit BTCS');

% ====== PLOT 3D ======
figure(2);
[T_mesh, X_mesh] = meshgrid(t, x);
surf(X_mesh, T_mesh, u, 'EdgeColor', 'none');
xlabel('x (m)');
ylabel('t (s)');
zlabel('u(x,t)');
title('Vibrasi Sasando 3D — Implisit BTCS');
colorbar;
view(30, 35);
