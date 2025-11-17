%ME4233 CA2
%Justin Kyle Leonard A0266108E
%Wang Zeqi A0306098Y

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

%Specify the value of dt
dt = 0.001; 
%The final time
tf = 2; 
%Maximum timestep
NMax = tf/dt;

%Assemble A matrix to solve the Poisson equation
dx2 = dx*dx; dy2 = dy*dy;
gamma2 = dx2/dy2;

[A] = Amatrix(Nx, Ny, dx2, gamma2);

%Initial condition for vorticity
vort = zeros(Nx-1, Ny-1);
t = 0;


%An array to record the kinetic energy
energyarray = [];

figure('units','normalized','outerposition',[0 0 1 1])

xc = 1.5;
yc = 0.5;

[~,ic] = min(abs(x-xc)); %Keep the index of the closest x-value to 1.5
[~,jc] = min(abs(y-yc)); %Keep the index of the closest y-value to 0.5

Uc = [];
Tc = [];

% time evoluion
for i=1:NMax

stmfunc = solvingpoisson(vort,A,Nx,Ny);           % solving the Poisson equation
                                                  % for streamfunction

vort = vortcalculation(stmfunc,vort,Nx,Ny,dx,dy,dt,Re,t); % advance the vorticity transport 
                                                          % equation                                                     

disp(['Finish the timestep ' num2str(i)])

t=t+dt;

% for plotting
figure(1)
subplot(1,2,1)
[u,v]= extract_uv(stmfunc,Nx,Ny,dx,dy,t);
uc_value = u(ic,jc);
Uc = [Uc; u(ic,jc)]; %Updating the list of uc
Tc = [Tc; t]; %Updating the list of lc

quiver(x2,y2,u',v')
axis equal
xlim([0 3]);ylim([0 1.1]);
title(['At time step = ' num2str(i)])
set(gca,'FontSize',40)
xlabel('x');ylabel('y')
pause(0.2)
% drawnow

energyarray = [energyarray sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];

subplot(1,2,2)
plot(dt:dt:dt*i,energyarray,'-*b')
title(['The kinetic energy at the grid point (' num2str(x(end-1)) ',' num2str(y(end-1)) ')'])
set(gca,'FontSize',20)
xlabel('time');ylabel('energy')

end

% post-processing
plot_stmfunc_vort(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2,t)

figure
subplot(1,2,1)
plot(Tc,Uc,'-b','LineWidth',1.5)    
title(['u velocity at point C(' num2str(xc) ',' num2str(yc) ') for Explicit Method'])    
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

%Implicit Method

%Initial condition for vorticity
vort = zeros(Nx-1, Ny-1);
ti = 0;
dti = 0.03;

NMaxi = tf/dti;

energyarrayim = []; %list of energy for implicit method

Uci = [];
Tci = [];

for i=1:NMaxi

stmfunc = solvingpoisson(vort,A,Nx,Ny);           % solving the Poisson equation
                                                  % for streamfunction
[LHS, RHS] = assembleLHSRHS(Nx,Ny,stmfunc,vort,dx,dy,dti,Re,ti);

LHS = sparse(LHS); RHS = sparse(RHS);
u0 = zeros((Nx-1)*(Ny-1),1);
resobj = 10^-4;
vort = SORmethod(LHS,RHS,u0,resobj); %Vorticity at the next timestep 
vort = reshape(vort,Nx-1,Ny-1);

disp(['Finish the timestep ' num2str(i)])

ti=ti+dti;

% for plotting
figure(5)
subplot(1,2,1)
[u,v]= extract_uv(stmfunc,Nx,Ny,dx,dy,ti);

Uci = [Uci; u(ic,jc)]; %Updating the list of uc
Tci = [Tci; ti]; %Updating the list of lc

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
plot_stmfunc_vort_implicit(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2,ti)

figure
subplot(1,2,1)
plot(Tci,Uci,'-b','LineWidth',1.5)    
title(['u velocity at point C(' num2str(xc) ',' num2str(yc) ') for Implicit Method'])    
set(gca,'FontSize',20)
xlabel('Time (t)')
ylabel('u velocity')

Fs = 1/dti;
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