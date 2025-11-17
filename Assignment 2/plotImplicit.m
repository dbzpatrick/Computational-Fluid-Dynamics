function [] = plotImplicit(stmfunc,vort,Nx,Ny,dx,dy,x,y,t)

% streamfunction
stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ];   
     
figure
contourf(x,y,stmfunc')
colorbar;
axis equal
set(gca,'FontSize',10)
title('Streamfunction at final timestep with Implicit Method')

% vorticity
dx2 = dx*dx; dy2 = dy*dy;
vortsouth=-2*stmfunc(2:end-1,2)/(dx2);
vortnorth=-2*stmfunc(2:end-1,end-1)/(dx2)-2*(sin(2*pi*5*t))/dx;
vorteast =-2*stmfunc(end-1,2:end-1)/(dy2);
vortwest =-2*stmfunc(2,2:end-1)/(dy2);
vort=[0         vortwest   0
      vortsouth vort       vortnorth
      0         vorteast   0 ];  
figure
contourf(x,y,vort')
colorbar;
axis equal  
set(gca,'FontSize',10)
title('Vorticity at final timestep with Implicit Method')
end