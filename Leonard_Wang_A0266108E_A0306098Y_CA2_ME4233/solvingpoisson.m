function [stmfunc] = solvingpoisson(vort, A, Nx, Ny)

b = -vort(:);

Ngs=(Nx-1)*(Ny-1);
u0=zeros(Ngs,1);
resobj=10^-4;

A=sparse(A);
stmfunc= SORmethod(A,b,u0,resobj);

% When we solve the Poisson equation, we stack the 2D grid points into a 1D
% vector. But when we solve the time-dependent equation, we can work with a
% 2D mesh. So here, we reshape the 1D vector back to the 2D mesh
stmfunc = reshape(stmfunc,Nx-1,Ny-1);

end