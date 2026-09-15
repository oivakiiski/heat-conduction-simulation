% Author: Oiva Kiiski
% Applied Mathematics 1 project work
% Topic: Heat conduction in solid medium
% Initial values:
% thermal conductivity λ = 205 W/m/K -> Aluminium
% Cylindrical segment: r = 0.4 m and l = 0.4 m
rhoo = 2701; % kg/m^3, density of aluminium at 300K
cp = 902; % J/kgK, specific heat capacity of aluminium at 300K
lambda = 205; % thermal conductivity in W/m/K
% thermal diffusivity is independent from geometry
alfa = lambda/(rhoo*cp); % m^2/s

R = 0.4; % radius in meters
L = 0.4; % length in meters
segments_r = 40;
segments_l = 40;
dr = R/segments_r; % spatial step in radial direction
dl = L/segments_l; % spatial step in axial direction

% nodes and initialization
Nr = segments_r + 1;
Nl = segments_l + 1;
r = linspace(0,R,Nr);
l = linspace(0,L,Nl);

% Time stepping and temperature
dt = 10;             % 1 second
tolerance = 0.000001;

T_hot = 60;
T_cold = 10;
T = T_cold*ones(Nr,Nl); % Initial temperature of 20 degrees
T(:,1) = T_hot;

% Build sparse laplacian operator for axisymmetric r-z domain
N = Nr*Nl;
index = @(i,j) (j-1)*Nr + i;
I = [];
J = [];
V = [];

% initialization of 5-point Laplacian

for j = 1:Nl
    for i = 1:Nr
        kindex = index(i,j);
        % Fixed boundary temperatures
        if (j==1) || (j==Nl) || (i==Nr)
            I(end+1) = kindex;
            J(end+1) = kindex;
            V(end+1) = 1;
            continue
        end

        ri = r(i);

        if i==1
            % The middle knot -> use ghost node
            L_diag = -2/dr^2 - 2/dl^2;
            L_ip = 2/dr^2;
            L_im = 0;
        else
            rp = ri + dr/2;
            rm = ri - dr/2;
            L_ip = (rp)/(ri*dr^2);
            L_im = (rm)/(ri*dr^2);
            L_diag = -(L_ip + L_im) - 2/dl^2;
        end
        
        % center
        I(end+1) = kindex;
        J(end+1) = kindex;
        V(end+1) = L_diag;
        % i+1
        I(end+1) = kindex;
        J(end+1) = index(i+1,j);
        V(end+1) = L_ip;
        % i-1
        if i>1
            I(end+1) = kindex;
            J(end+1) = index(i-1,j);
            V(end+1) = L_im;
        end
        % j+1
        I(end+1)=kindex;
        J(end+1) = index(i,j+1);
        V(end+1) = 1/dl^2;
        % j-1
        I(end+1)=kindex;
        J(end+1) = index(i,j-1);
        V(end+1) = 1/dl^2;
    end
end

Laplacian = sparse(I,J,V,N,N);
% Build CN system matrices
A = speye(N) - 1/2*alfa*dt*Laplacian;
B = speye(N) + 1/2*alfa*dt*Laplacian;

% Physical time step
Tvector = T(:);
t = 0;
step = 0;
% Preset boundary temperatures
BC_index = false(N,1);
BC_values = zeros(N,1);
for j = 1:Nl
    for i = 1:Nr
        kindex = index(i,j);
        if j==1
            BC_index(kindex) = true;
            BC_values(kindex) = T_hot; % Set boundary condition for hot side
        elseif (j==Nl) || (i==Nr)
            BC_index(kindex) = true;
            BC_values(kindex) = T_cold; % Set boundary condition for cold side
        end
    end
end

% Crank-Nicolson: ATn+1=BTn

while true
    step = step + 1;
    t = t + dt;
    T_old = Tvector;
    rhs = B*Tvector;
    % Boundary conditions into rhs
    rhs(BC_index) = BC_values(BC_index);
    % Solve linear system A*Tn+1=rhs
    Tvector = A \ rhs;
    % tolerance check
    if max(abs(Tvector-T_old)) < tolerance
        break
    end
end
fprintf('Steady state reched at t = %.2f s (steps = %d)\n',t,step);

% Resaphe into 2D field
T = reshape(Tvector, Nr, Nl);
% Central difference & Fourier's law. Heat flux
qr = zeros(Nr, Nl);
ql = zeros(Nr, Nl);
for j = 2:Nl-1
    for i = 2:Nr-1
        dTdr = (T(i+1,j) - T(i-1,j)) / (2*dr);
        dTdl = (T(i,j+1) - T(i,j-1)) / (2*dl);
        qr(i,j) = -lambda * dTdr;
        ql(i,j) = -lambda * dTdl;
    end
    % axis i=1
    dTdr_axis = (T(2,j) - T(1,j)) / dr;
    qr(1,j) = -lambda * dTdr_axis;
end

% plots
[Z, Rg] = meshgrid(r, l);
figure; 
contourf(Rg, Z, T, 30, 'EdgeColor','none'); 
colorbar;
xlabel('r (m)'); 
ylabel('l (m)');
title('Steady-State Temperature (°C) — Crank–Nicolson');

% quiver plot of heat flux
skip = 3;
figure; 
quiver(Rg(1:skip:end,1:skip:end)', Z(1:skip:end,1:skip:end)', qr(1:skip:end,1:skip:end)', ql(1:skip:end,1:skip:end)');
xlabel('r (m)'); 
ylabel('l (m)'); 
title('Heat flux vectors (q_l, q_r)');
