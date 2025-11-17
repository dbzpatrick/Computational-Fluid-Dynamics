%Implicit Method

clear all;

%Specify Reynolds Number
Re = 28; 

%Set up the grids
Nx = 50; Ny = 20;
Lx = 3; Ly = 1;
x = linspace(0, Lx, Nx+1);
y = linspace(0, Ly, Ny+1);
dx = x(2) - x(1);
dy = y(2) - y(1);

[x2, y2] = meshgrid(x,y);

%Assemble A matrix to solve the Poisson equation
dx2 = dx*dx; dy2 = dy*dy;
gamma2 = dx2/dy2;

[A] = Amatrix(Nx, Ny, dx2, gamma2);


%Initial condition for vorticity
vorti = zeros(Nx-1, Ny-1);
ti = 0;
tf = 2;
dti = 0.05;

NMaxi = tf/dti;

energyarrayim = []; %list of energy uhhuyhuifor implicit method


for i=1:NMaxi

stmfunc = solvingpoisson(vorti,A,Nx,Ny);           % solving the Poisson equation
                                                  % for streamfunction
vorti = vortcalculation_implicit(Nx,Ny,stmfunc,vorti,dx,dy,dti,Re,ti);

disp(['Finish the timestep ' num2str(i)])

ti=ti+dti;

% for plotting
figure(1)
subplot(1,2,1)
[u,v]= extract_uv(stmfunc,Nx,Ny,dx,dy,ti);
quiver(x2,y2,u',v')
axis equal
xlim([0 3]);ylim([0 1.1]);
title(['At time step = ' num2str(i)])
set(gca,'FontSize',40)
xlabel('x');ylabel('y')
pause(0.2)
% drawnow

energyarrayim = [energyarrayim sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];

subplot(1,2,2)
plot(dti:dti:dti*i,energyarrayim,'-*b')
title(['The kinetic energy at the grid point (' num2str(x(end-1)) ',' num2str(y(end-1)) ')'])
set(gca,'FontSize',20)
xlabel('time');ylabel('energy')

end

% post-processing
plot_stmfunc_vort_implicit(stmfunc,vorti,Nx,Ny,dx,dy,x,y,x2,y2,ti)

figure
subplot(1,2,1)
plot(Tc,Uc,'-b','LineWidth',1.5)    
title(['u velocity at point C(' num2str(xc) ',' num2str(yc) ') for Implicit Method'])    
set(gca,'FontSize',20)
xlabel('Time (t)')
ylabel('u velocity')

Fs = 1/dt;
L = length(Uc);
Y = fft(Uc);
P2 = abs(Y/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);
ff = Fs*(0:(L-1)/2)/L;

subplot(1,2,2)
plot(ff,P1,"-o") 
title("Single-Sided Amplitude")
xlabel("f (Hz)")
ylabel("|P1(f)|")
set(gca,'FontSize',20)
