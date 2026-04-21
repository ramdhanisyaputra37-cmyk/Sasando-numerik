clc; clear all; close all;

% ====== PARAMETER ======
c  = 1;       % kecepatan gelombang
kd = 1.5;     % koefisien redaman
l  = 0.64;    % panjang senar
h  = 1;       % amplitudo petikan
dx = 0.05;    % step spasial
dt = 0.008;     % step waktu

% ====== GRID ======
x = 0:dx:l;
t = 0:dt:7;
Nx = length(x);
Nt = length(t);

% ====== KONDISI AWAL (petikan segitiga) ======
u = zeros(Nx, Nt);
for j = 1:Nx
    if x(j) <= 0.5*l
        u(j,1) = 2*h*x(j)/l;
    else
        u(j,1) = 2*h - 2*h*x(j)/l;
    end
end

% ====== PLOT KONDISI AWAL ======
figure(1);
plot(x, u(:,1), 'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('u(x,0)');
title('Kondisi Awal — Petikan Senar Sasando');
grid on;
% ====== KONDISI BATAS ======
% Ujung kiri dan kanan senar selalu = 0
for n = 1:Nt
    u(1, n)  = 0;   % u(0,t) = 0
    u(Nx, n) = 0;   % u(l,t) = 0
end

% ====== CEK KONDISI BATAS ======
fprintf('Ujung kiri  u(1,1)  = %.4f\n', u(1,1));
fprintf('Ujung kanan u(Nx,1) = %.4f\n', u(Nx,1));
% ====== KOEFISIEN A ======
A = 0.5*c^2 + 2*c^2/l;
r = A * dt^2 / dx^2;

fprintf('A = %.4f\n', A);
fprintf('r = %.4f\n', r);

% ====== INISIALISASI LEVEL n=1 ======
% Dari du/dt(x,0) = 0  →  u^{-1} = u^0
% Sehingga: u(j,2) = u(j,1) + 0.5*r*(u(j+1,1) - 2*u(j,1) + u(j-1,1))
for j = 2:Nx-1
    u(j,2) = u(j,1) + 0.5*r*(u(j+1,1) - 2*u(j,1) + u(j-1,1));
end

% ====== PLOT PERBANDINGAN n=0 dan n=1 ======
figure(2);
plot(x, u(:,1), 'b-', 'LineWidth', 2); hold on;
plot(x, u(:,2), 'r--', 'LineWidth', 2);
xlabel('x (m)');
ylabel('u(x,t)');
title('Kondisi Awal vs Level n=1');
legend('n=0 (t=0)', 'n=1 (t=dt)');
grid on;
% ====== ITERASI EKSPLISIT KLASIK ======
D = 1 + kd*dt/2;

for n = 2:Nt-1
    for j = 2:Nx-1
        lap = u(j+1,n) - 2*u(j,n) + u(j-1,n);
        vel = u(j,n) - u(j,n-1);
        u(j,n+1) = (2*u(j,n) - u(j,n-1) + r*lap - kd*dt*vel) / D;
    end
end

% ====== PLOT HASIL ======
figure(3);
for n = 1:Nt
    plot(x, u(:,n), 'b-', 'LineWidth', 1.5);
    xlabel('x (m)');
    ylabel('u(x,t)');
    title('Vibrasi Sasando — Eksplisit Klasik');
    grid on;
    hold on;
end
% ====== PLOT 3D ======
figure(4);
[T_mesh, X_mesh] = meshgrid(t, x);
surf(X_mesh, T_mesh, u, 'EdgeColor', 'none');
xlabel('x (m)');
ylabel('t (s)');
zlabel('u(x,t)');
title('Vibrasi Sasando 3D — Eksplisit Klasik');
colorbar;
view(30, 35);
