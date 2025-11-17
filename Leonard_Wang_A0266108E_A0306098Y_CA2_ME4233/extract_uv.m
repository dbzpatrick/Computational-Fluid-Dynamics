function [u,v] = extract_uv(stmfunc,Nx,Ny,dx,dy,t)
%GET_UV Summary of this function goes here
%   Detailed explanation goes here

% Surround the streamfunction (interior grid points) with the bc
stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ]; 
     
% Get u,v from the streamfunction
u =  (stmfunc(2:end-1,3:end) - stmfunc(2:end-1,1:end-2))/(2*dx);
v = -(stmfunc(3:end,2:end-1) - stmfunc(1:end-2,2:end-1))/(2*dy);

% Surround the u,v (interior grid points) with the bc
u=[0             zeros(1,Ny-1)  sin(2*pi*5*t)
   zeros(Nx-1,1) u              ones(Nx-1,1)*sin(2*pi*5*t)
   0             zeros(1,Ny-1)  sin(2*pi*5*t) ]; 

v=[0             zeros(1,Ny-1)  0
   zeros(Nx-1,1) v              zeros(Nx-1,1)
   0             zeros(1,Ny-1)  0 ]; 

end

