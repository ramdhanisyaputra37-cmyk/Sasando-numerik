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

% ====== KOEFISIEN ======
A   = 0.5*c^2 + 2*c^2/l;
r   = A*dt^2/dx^2;
cw  = sqrt(A);
mu  = A*dt/dx^2;

% ====== KONDISI AWAL ======
u_eks = zeros(Nx,Nt);
u_imp = zeros(Nx,Nt);
u_df  = zeros(Nx,Nt);

for j = 1:Nx
    if x(j) <= 0.5*l
        val = 2*h*x(j)/l;
    else
        val = 2*h - 2*h*x(j)/l;
    end
    u_eks(j,1) = val;
    u_imp(j,1) = val;
    u_df(j,1)  = val;
end

% ====== KONDISI BATAS ======
for n = 1:Nt
    u_eks(1,n)=0; u_eks(Nx,n)=0;
    u_imp(1,n)=0; u_imp(Nx,n)=0;
    u_df(1,n)=0;  u_df(Nx,n)=0;
end

% ====== INISIALISASI n=1 ======
for j = 2:Nx-1
    u_eks(j,2) = u_eks(j,1) + 0.5*r*(u_eks(j+1,1)-2*u_eks(j,1)+u_eks(j-1,1));
    u_imp(j,2) = u_imp(j,1) + 0.5*r*(u_imp(j+1,1)-2*u_imp(j,1)+u_imp(j-1,1));
    u_df(j,2)  = u_df(j,1)  + 0.5*r*(u_df(j+1,1) -2*u_df(j,1) +u_df(j-1,1));
end

% ====== ITERASI EKSPLISIT ======
D_eks = 1 + kd*dt/2;
for n = 2:Nt-1
    for j = 2:Nx-1
        lap = u_eks(j+1,n)-2*u_eks(j,n)+u_eks(j-1,n);
        vel = u_eks(j,n)-u_eks(j,n-1);
        u_eks(j,n+1) = (2*u_eks(j,n)-u_eks(j,n-1)+r*lap-kd*dt*vel)/D_eks;
    end
end

% ====== ITERASI IMPLISIT ======
Da = 2*l*dx^2+2*cw^2*l*dt^2+8*cw^2*dt^2+2*kd*l*dx^2*dt;
Db = l*dx^2+cw^2*l*dt^2+4*cw^2*dt^2+kd*l*dx^2*dt;
aa = (cw^2*l*dt^2+4*cw^2*dt^2)/Da;
bb = (2*l*dx^2+kd*l*dx^2*dt)/Db;
cc = l*dx^2/Db;
M     = Nx-2;
A_tri = diag(ones(M,1))-aa*diag(ones(M-1,1),1)-aa*diag(ones(M-1,1),-1);
for n = 2:Nt-1
    B = bb*u_imp(2:Nx-1,n) - cc*u_imp(2:Nx-1,n-1);
    u_imp(2:Nx-1,n+1) = A_tri \ B;
end

% ====== ITERASI DUFORT-FRANKLE ======
D_df = 1+2*mu+2*kd*dt;
for n = 2:Nt-1
    for j = 2:Nx-1
        neigh = u_df(j+1,n)+u_df(j-1,n);
        u_df(j,n+1) = ((1-2*mu)*u_df(j,n-1)+2*mu*neigh-2*kd*dt*u_df(j,n))/D_df;
    end
end

% ====== PLOT PERBANDINGAN SNAPSHOT ======
t_snap = [1, 40, 80, 120, 160, 188];   % indeks waktu
t_label = {'t=0.0', 't=0.3', 't=0.6', 't=0.9', 't=1.2', 't=1.5'};

figure(1);
for k = 1:6
    subplot(2,3,k);
    n = t_snap(k);
    plot(x, u_eks(:,n), 'r-',  'LineWidth', 2); hold on;
    plot(x, u_imp(:,n), 'b--', 'LineWidth', 2);
    plot(x, u_df(:,n),  'g:',  'LineWidth', 2);
    xlabel('x (m)'); ylabel('u(x,t)');
    title(t_label{k});
    ylim([-1.3 1.3]);
    grid on;
    if k == 1
        legend('Eksplisit','Implisit','Dufort-F.','Location','northeast');
    end
end


% ====== PLOT EVOLUSI TITIK TENGAH ======
mid = round(Nx/2);
pause(0.1);
figure(2);
plot(t, u_eks(mid,:), 'r-',  'LineWidth', 2); hold on;
plot(t, u_imp(mid,:), 'b--', 'LineWidth', 2);
plot(t, u_df(mid,:),  'g:',  'LineWidth', 2);
xlabel('t (s)'); ylabel('u');
title('Evolusi Waktu di Titik Tengah x=0.5l');
legend('Eksplisit','Implisit','Dufort-Frankle');
grid on;
% ====== PLOT 3D PERBANDINGAN ======
[T_mesh, X_mesh] = meshgrid(t, x);

figure(3);
subplot(1,3,1);
surf(X_mesh, T_mesh, u_eks, 'EdgeColor','none');
xlabel('x (m)'); ylabel('t (s)'); zlabel('u');
title('Eksplisit Klasik');
colorbar; view(30,35);

subplot(1,3,2);
surf(X_mesh, T_mesh, u_imp, 'EdgeColor','none');
xlabel('x (m)'); ylabel('t (s)'); zlabel('u');
title('Implisit BTCS');
colorbar; view(30,35);

subplot(1,3,3);
surf(X_mesh, T_mesh, u_df, 'EdgeColor','none');
xlabel('x (m)'); ylabel('t (s)'); zlabel('u');
title('Dufort-Frankle');
colorbar; view(30,35);
% ====== PLOT HEATMAP ======
figure(4);

% Baris 1: Solusi masing-masing metode
subplot(2,3,1);
pcolor(T_mesh, X_mesh, u_eks); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('Eksplisit Klasik');
colorbar;

subplot(2,3,2);
pcolor(T_mesh, X_mesh, u_imp); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('Implisit BTCS');
colorbar;

subplot(2,3,3);
pcolor(T_mesh, X_mesh, u_df); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('Dufort-Frankle');
colorbar;

% Baris 2: Selisih antar metode
subplot(2,3,4);
pcolor(T_mesh, X_mesh, abs(u_imp - u_eks)); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('|Implisit - Eksplisit|');
colorbar;

subplot(2,3,5);
pcolor(T_mesh, X_mesh, abs(u_imp - u_df)); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('|Implisit - Dufort-F.|');
colorbar;

subplot(2,3,6);
pcolor(T_mesh, X_mesh, abs(u_eks - u_df)); shading interp;
xlabel('t (s)'); ylabel('x (m)');
title('|Eksplisit - Dufort-F.|');
colorbar;
