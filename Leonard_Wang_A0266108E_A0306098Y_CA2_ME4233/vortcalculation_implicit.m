function [vortnew] = vortcalculation_implicit(stmfunc,vorti,Nx,Ny,dx, dy ,dti,Re, ti)

[LHS, RHS] = assembleLHSRHS(Nx, Ny, stmfunc, vorti, dx, dy ,dti,Re, ti);

LHS = sparse(LHS); RHS = sparse(RHS);
Ngs = (Nx-1)*(Ny-1);
u0 = zeros(Ngs,1);
resobj = 10^-4;

vortnew = SORmethod(LHS,RHS,Ngs, u0, resobj);