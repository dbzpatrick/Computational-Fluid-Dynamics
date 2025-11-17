% function [RHS] = assembleRHS(Nx,Ny,stmfunc,vort,Re,dx,dy)
% %ASSEMBLERHS Summary of this function goes here
% %   Detailed explanation goes here
% 
% % the boundary conditions at the boundaries
% U_south = zeros(Nx-1,1);
% U_north = ones(Nx-1,1);   % This is the driving force
% U_west  = zeros(Ny-1,1);
% U_east  = zeros(Ny-1,1);
% 
% nu=1/Re;
% h1=dx*dx;
% h2=dy*dy;
% 
% RHS=zeros(Nx-1,Ny-1);
% 
% % interior
% for j=2:Ny-2
%     for i=2:Nx-2
%         fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
%         fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
% 
%         RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
%                     fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
%                     nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                          ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/h2 );                     
%     end
% end
% 
% % south
% j=1;
% 
% for i=2:Nx-2
%     vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*h2);
% 
%     fac1 = -(stmfunc(i,j+1) - 0             )/2/dy;
%     fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
% 
%     RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
%                 fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
%                 nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                      ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/h2 ); 
% end
% 
% % North
% j=Ny-1;
% for i=2:Nx-2
%     vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*h2);
% 
%     fac1 = -(0              - stmfunc(i,j-1))/2/dy;
%     fac2 =  (stmfunc(i+1,j) - stmfunc(i-1,j))/2/dx;
% 
%     RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
%                 fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
%                 nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                      ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/h2 ); 
% end
% 
% % West
% i=1;
% for j=2:Ny-2
%     vortwestbc = ( 0 - stmfunc(i,j) - U_west(j)*dx)/(0.5*h1);
% 
%     fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
%     fac2 =  (stmfunc(i+1,j) - 0             )/2/dx;
% 
%     RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
%                 fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
%                 nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/h1 + ... 
%                      ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/h2 ); 
% end
% 
% % East
% i=Nx-1;
% for j=2:Ny-2
%     vorteastbc = ( 0 - stmfunc(i,j) + U_east(j)*dx)/(0.5*h1);
% 
%     fac1 = -(stmfunc(i,j+1) - stmfunc(i,j-1))/2/dy;
%     fac2 =  (0              - stmfunc(i-1,j))/2/dx;
% 
%     RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
%                 fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
%                 nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                      ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/h2 ); 
% end
% 
% % South-west
% i=1;j=1;
% 
%     vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*h2);
%     vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*h1);
% 
%     fac1 = -(stmfunc(i,j+1) - 0 )/2/dy;
%     fac2 =  (stmfunc(i+1,j) - 0 )/2/dx;
% 
%     RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
%                 fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
%                 nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/h1 + ... 
%                      ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/h2 ); 
% 
% % South-east          
% i=Nx-1;j=1;
% 
%     vortsouthbc = ( 0 - stmfunc(i,j) + U_south(i)*dy)/(0.5*h2);
%     vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*h1);
% 
%     fac1 = -(stmfunc(i,j+1) - 0)/2/dy;
%     fac2 =  (0 - stmfunc(i-1,j))/2/dx;
% 
%     RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
%                 fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
%                 nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                      ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/h2 ); 
% 
% % North-east         
% i=Nx-1;j=Ny-1;
% 
%     vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*h2);
%     vorteastbc  = ( 0 - stmfunc(i,j) - U_east(j) *dx)/(0.5*h1);
% 
%     fac1 = -(0 - stmfunc(i,j-1))/2/dy;
%     fac2 =  (0 - stmfunc(i-1,j))/2/dx;
% 
%     RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
%                 fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
%                 nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/h1 + ... 
%                      ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/h2 );
% 
% % North-west          
% i=1;j=Ny-1;
% 
%     vortnorthbc = ( 0 - stmfunc(i,j) - U_north(i)*dy)/(0.5*h2);
%     vortwestbc  = ( 0 - stmfunc(i,j) + U_west(j) *dx)/(0.5*h1);
% 
%     fac1 = -(0 - stmfunc(i,j-1))/2/dy;
%     fac2 =  (stmfunc(i+1,j) - 0)/2/dx;
% 
%     RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
%                 fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
%                 nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/h1 + ... 
%                      ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/h2 ); 
% 
% end
% 
% function vortnew = advance_vort(stmfunc, vort, Nx, Ny, dx, dy, dt, Re)
% % Advances the vorticity field by one time step
% % Uses explicit Euler for time integration
% %
% % Inputs:
% %   stmfunc: (Ny+1) x (Nx+1) - streamfunction INCLUDING boundaries
% %   vort:    (Ny-1) x (Nx-1) - vorticity at INTERIOR points only
% %   Nx, Ny:  number of intervals
% %   dx, dy:  grid spacing
% %   dt:      time step
% %   Re:      Reynolds number
% 
% RHS = assembleRHS(Nx, Ny, stmfunc, vort, Re, dx, dy);
% 
% % Explicit Euler method
% vortnew = vort + dt * RHS;
% 
% end

function RHS = assembleRHS(Nx, Ny, stmfunc, vort, Re, dx, dy)
% Your original code structure but with fixes
% Matches YOUR indexing convention: 
%   vort(i,j) and stmfunc(i,j) where i=1:Nx-1, j=1:Ny-1

nu = 1/Re;
h1 = dx^2;
h2 = dy^2;

% Boundary velocities
U_south = zeros(Nx-1, 1);
U_north = ones(Nx-1, 1);   % Moving lid
U_west  = zeros(Ny-1, 1);
U_east  = zeros(Ny-1, 1);

RHS = zeros(Nx-1, Ny-1);

%% Interior points
for j = 2:Ny-2
    for i = 2:Nx-2
        % Velocities: u = ∂ψ/∂y, v = -∂ψ/∂x
        fac1 = -(stmfunc(i, j+1) - stmfunc(i, j-1)) / (2*dy);  % u
        fac2 =  (stmfunc(i+1, j) - stmfunc(i-1, j)) / (2*dx);  % v
        
        RHS(i,j) = fac1 * (vort(i+1,j) - vort(i-1,j)) / (2*dx) + ...
                   fac2 * (vort(i,j+1) - vort(i,j-1)) / (2*dy) + ...
                   nu * ((vort(i+1,j) - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                         (vort(i,j+1) - 2*vort(i,j) + vort(i,j-1)) / h2);
    end
end

%% South boundary (j=1)
j = 1;
for i = 2:Nx-2
    % Boundary vorticity using no-slip condition
    vortsouthbc = -2*stmfunc(i,j) / dy^2 + 2*U_south(i) / dy;
    
    fac1 = -(stmfunc(i, j+1) - 0) / (2*dy);
    fac2 =  (stmfunc(i+1, j) - stmfunc(i-1, j)) / (2*dx);
    
    RHS(i,j) = fac1 * (vort(i+1,j) - vort(i-1,j)) / (2*dx) + ...
               fac2 * (vort(i,j+1) - vortsouthbc) / (2*dy) + ...
               nu * ((vort(i+1,j) - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                     (vort(i,j+1) - 2*vort(i,j) + vortsouthbc) / h2);
end

%% North boundary (j=Ny-1)
j = Ny-1;
for i = 2:Nx-2
    vortnorthbc = -2*stmfunc(i,j) / dy^2 - 2*U_north(i) / dy;
    
    fac1 = -(0 - stmfunc(i, j-1)) / (2*dy);
    fac2 =  (stmfunc(i+1, j) - stmfunc(i-1, j)) / (2*dx);
    
    RHS(i,j) = fac1 * (vort(i+1,j) - vort(i-1,j)) / (2*dx) + ...
               fac2 * (vortnorthbc - vort(i,j-1)) / (2*dy) + ...
               nu * ((vort(i+1,j) - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                     (vortnorthbc - 2*vort(i,j) + vort(i,j-1)) / h2);
end

%% West boundary (i=1)
i = 1;
for j = 2:Ny-2
    vortwestbc = -2*stmfunc(i,j) / dx^2 + 2*U_west(j) / dx;
    
    fac1 = -(stmfunc(i, j+1) - stmfunc(i, j-1)) / (2*dy);
    fac2 =  (stmfunc(i+1, j) - 0) / (2*dx);
    
    RHS(i,j) = fac1 * (vort(i+1,j) - vortwestbc) / (2*dx) + ...
               fac2 * (vort(i,j+1) - vort(i,j-1)) / (2*dy) + ...
               nu * ((vort(i+1,j) - 2*vort(i,j) + vortwestbc) / h1 + ...
                     (vort(i,j+1) - 2*vort(i,j) + vort(i,j-1)) / h2);
end

%% East boundary (i=Nx-1)
i = Nx-1;
for j = 2:Ny-2
    vorteastbc = -2*stmfunc(i,j) / dx^2 - 2*U_east(j) / dx;
    
    fac1 = -(stmfunc(i, j+1) - stmfunc(i, j-1)) / (2*dy);
    fac2 =  (0 - stmfunc(i-1, j)) / (2*dx);
    
    RHS(i,j) = fac1 * (vorteastbc - vort(i-1,j)) / (2*dx) + ...
               fac2 * (vort(i,j+1) - vort(i,j-1)) / (2*dy) + ...
               nu * ((vorteastbc - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                     (vort(i,j+1) - 2*vort(i,j) + vort(i,j-1)) / h2);
end

%% Corner: South-west
i = 1; j = 1;
vortsouthbc = -2*stmfunc(i,j) / dy^2 + 2*U_south(i) / dy;
vortwestbc  = -2*stmfunc(i,j) / dx^2 + 2*U_west(j) / dx;

fac1 = -(stmfunc(i, j+1) - 0) / (2*dy);
fac2 =  (stmfunc(i+1, j) - 0) / (2*dx);

RHS(i,j) = fac1 * (vort(i+1,j) - vortwestbc) / (2*dx) + ...
           fac2 * (vort(i,j+1) - vortsouthbc) / (2*dy) + ...
           nu * ((vort(i+1,j) - 2*vort(i,j) + vortwestbc) / h1 + ...
                 (vort(i,j+1) - 2*vort(i,j) + vortsouthbc) / h2);

%% Corner: South-east
i = Nx-1; j = 1;
vortsouthbc = -2*stmfunc(i,j) / dy^2 + 2*U_south(i) / dy;
vorteastbc  = -2*stmfunc(i,j) / dx^2 - 2*U_east(j) / dx;

fac1 = -(stmfunc(i, j+1) - 0) / (2*dy);
fac2 =  (0 - stmfunc(i-1, j)) / (2*dx);

RHS(i,j) = fac1 * (vorteastbc - vort(i-1,j)) / (2*dx) + ...
           fac2 * (vort(i,j+1) - vortsouthbc) / (2*dy) + ...
           nu * ((vorteastbc - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                 (vort(i,j+1) - 2*vort(i,j) + vortsouthbc) / h2);

%% Corner: North-east
i = Nx-1; j = Ny-1;
vortnorthbc = -2*stmfunc(i,j) / dy^2 - 2*U_north(i) / dy;
vorteastbc  = -2*stmfunc(i,j) / dx^2 - 2*U_east(j) / dx;

fac1 = -(0 - stmfunc(i, j-1)) / (2*dy);
fac2 =  (0 - stmfunc(i-1, j)) / (2*dx);

RHS(i,j) = fac1 * (vorteastbc - vort(i-1,j)) / (2*dx) + ...
           fac2 * (vortnorthbc - vort(i,j-1)) / (2*dy) + ...
           nu * ((vorteastbc - 2*vort(i,j) + vort(i-1,j)) / h1 + ...
                 (vortnorthbc - 2*vort(i,j) + vort(i,j-1)) / h2);

%% Corner: North-west
i = 1; j = Ny-1;
vortnorthbc = -2*stmfunc(i,j) / dy^2 - 2*U_north(i) / dy;
vortwestbc  = -2*stmfunc(i,j) / dx^2 + 2*U_west(j) / dx;

fac1 = -(0 - stmfunc(i, j-1)) / (2*dy);
fac2 =  (stmfunc(i+1, j) - 0) / (2*dx);

RHS(i,j) = fac1 * (vort(i+1,j) - vortwestbc) / (2*dx) + ...
           fac2 * (vortnorthbc - vort(i,j-1)) / (2*dy) + ...
           nu * ((vort(i+1,j) - 2*vort(i,j) + vortwestbc) / h1 + ...
                 (vortnorthbc - 2*vort(i,j) + vort(i,j-1)) / h2);

end