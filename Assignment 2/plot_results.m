function [] = plot_results(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2)

% plot the streamfunction with boundary conditions
stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ];        
figure
contourf(x,y,stmfunc')
axis equal
set(gca,'FontSize',10)
title('Streamfunction at final time step with Explicit Method')

% plot the vorticity  with boundary conditions
vortsouth=-2*stmfunc(2:end-1,2)/dy^2;
vortnorth=-2*stmfunc(2:end-1,end-1)/dy^2-2*1/dy; % 1 in the last term means U_north=1
vorteast =-2*stmfunc(end-1,2:end-1)/dx^2;
vortwest =-2*stmfunc(2,2:end-1)/dx^2;

vort=[0         vortwest   0
      vortsouth vort       vortnorth
      0         vorteast   0 ];  
figure
contourf(x,y,vort')
axis equal  
set(gca,'FontSize',10)
title('Vorticity at final time step with Explicit Method')
end

