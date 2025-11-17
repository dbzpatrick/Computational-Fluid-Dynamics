function [] = plot_stmfunc_vort_implicit(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2, t)
%PLOT_RESULTS Summary of this function goes here
%   Detailed explanation goes here

% retrieve the velocity field from the streamfunction
[u,v]=extract_uv(stmfunc,Nx,Ny,dx,dy,t);

%Reshape the stmfunc and vort

% plot the streamfunction with boundary conditions
stmfunc=[0             zeros(1,Ny-1)   0
         zeros(Nx-1,1) stmfunc         zeros(Nx-1,1)
         0             zeros(1,Ny-1)   0 ];   
     
figure
contourf(x,y,stmfunc')
colorbar;
axis equal
set(gca,'FontSize',40)
title('Streamfunction at the final timestep (Implicit)')

% plot the vorticity  with boundary conditions
dx2 = dx*dx; dy2 = dy*dy;
vortsouth=-2*stmfunc(2:end-1,2)/(dx2);
vortnorth=-2*stmfunc(2:end-1,end-1)/(dx2)-2*(sin(2*pi*5*t))/dx; % 1 in the last term means U_north=1
vorteast =-2*stmfunc(end-1,2:end-1)/(dy2);
vortwest =-2*stmfunc(2,2:end-1)/(dy2);

vort=[0         vortwest   0
      vortsouth vort       vortnorth
      0         vorteast   0 ];  

figure
contourf(x,y,vort')
colorbar;
axis equal  
set(gca,'FontSize',40)
title('Vorticity at the final timestep (Implicit)')

% % plot the velocity field with the boundary conditions
% u = (u(:,2:end)+u(:,1:end-1))/2;
% u = (u(2:end,:)+u(1:end-1,:))/2;
% 
% v = (v(:,2:end)+v(:,1:end-1))/2;
% v = (v(2:end,:)+v(1:end-1,:))/2;
% 
% x2=(x2(2:end,2:end)+x2(2:end,1:end-1))/2;
% y2=(y2(2:end,2:end)+y2(1:end-1,2:end))/2;
% 
% figure
% quiver(x2,y2,u',v')
% axis equal
% set(gca,'FontSize',40)


end