function [vortnew] = vortcalculation(stmfunc,vort,Nx,Ny,dx, dy ,dt,Re, t)
%SOLVE_ Summary of this function goes here
%   Detailed explanation goes here

RHS = RHSassembly(Nx, Ny, stmfunc, vort, Re, dx, dy, t);

% Explicit Euler method
vortnew = vort + dt*RHS;
                                
end