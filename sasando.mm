clc; clear all; close all;

% ====== PARAMETER (sama dengan PPT) ======
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

% ====== KOEFISIEN ======
A = 0.5*c^2 + 2*c^2/l;
r = A*dt^2/dx^2;
fprintf('A = %.4f\n', A);
fprintf('r = %.4f\n', r);

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

% ====== ITERASI EKSPLISIT KLASIK ======
D = 1 + kd*dt/2;
for n = 2:Nt-1
    for j = 2:Nx-1
        lap = u(j+1,n) - 2*u(j,n) + u(j-1,n);
        vel = u(j,n) - u(j,n-1);
        u(j,n+1) = (2*u(j,n) - u(j,n-1) + r*lap - kd*dt*vel) / D;
    end
end

% ====== PLOT 2D ======
figure(1);
hold on; grid on;
for n = 1:Nt
    plot(x, u(:,n), 'b-', 'LineWidth', 1);
end
xlabel('x (m)');
ylabel('u(x,t)');
title('Vibrasi Sasando — Eksplisit Klasik');

% ====== PLOT 3D ======
figure(2);
[T_mesh, X_mesh] = meshgrid(t, x);
surf(X_mesh, T_mesh, u, 'EdgeColor', 'none');
xlabel('x (m)');
ylabel('t (s)');
zlabel('u(x,t)');
title('Vibrasi Sasando 3D — Eksplisit Klasik');
colorbar;
view(30, 35);
