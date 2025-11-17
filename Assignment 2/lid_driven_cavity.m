%% ME4233 Part 2 Assignment 2

% Dai Baizhou A0266515B
% Sam Marret A0252323R
% Aditiya Satish Nalini A0244516H


%% part A

%  A lid-driven cavity flow problem using Streamfunction-Vorticity formulation
clear;

% specify Reynolds number
Re=25; 

% set up the grids
Nx=51;Ny=21;
Lx=3;Ly=1;
dx=Lx/Nx; dy=Ly/Ny;
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
[x2,y2]=meshgrid(x,y);  

% specify the value of dt
dt=0.001;   
t_f=0.02;
NMax=t_f/dt;

% Assemble the A matrix to solve the Poisson equation Au=b
gamma = dx^2/dy^2;
A=Amatrix(Nx,Ny,dx^2,gamma);

% initial condition for the vorticity
vort = zeros(Nx-1,Ny-1);
t=0;

% point C
Cx = 1.5; Cy = 0.5;
[~,ic] = min(abs(x-Cx));
[~,jc] = min(abs(y-Cy));
uC = [];
timeC = [];

% An array to record the kinetic energy at a certain point
energyarray=[];

% time evoluion
for i=1:NMax
stmfunc = solve_Poisson(vort,A,Nx,Ny);            % solving the Poisson equation
vort = advance_vort(stmfunc,vort,Nx,Ny,dx,dy,dt,Re); % advance the vorticity transport 
disp(['Finish the timestep ' num2str(i)])
t=t+dt;

% for plotting
figure(1)
subplot(1,2,1)
[u,v]=get_uv(stmfunc,Nx,Ny,dx,dy);
uC = [uC; u(ic,jc)];timeC = [timeC,t];
quiver(x2,y2,u',v',5)
title(['At time step = ' num2str(i)])
set(gca,'FontSize',10)
xlabel('x');ylabel('y');ylim([0,1]);xlim([0,3]);
pause(0.2)
% drawnow
energyarray = [energyarray sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];
subplot(1,2,2)
plot(dt:dt:dt*i,energyarray,'-*b')
title(['The kinetic energy at the grid point (' num2str(x(end-1)) ',' num2str(y(end-1)) ')'])
set(gca,'FontSize',10)
xlabel('time');ylabel('energy')
end

%% post-processing
plot_results(stmfunc,vort,Nx,Ny,dx,dy,x,y,x2,y2)
figure
subplot(1,2,1)
plot(timeC,uC,'-b','LineWidth',1.5)    
title(['u velocity at point C(1.5,0.5) for Explicit Method'])    
set(gca,'FontSize',10)
xlabel('time')
ylabel('u signal')
Fs = 1/dt;
L = length(uC);
Y = fft(uC);
P2 = abs(Y/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);
ff = Fs*(0:(L-1)/2)/L;
subplot(1,2,2)
plot(ff,P1,"-x") 
title("Single-Sided Amplitude")
xlabel("f (Hz)")
ylabel("P1 (f)")
set(gca,'FontSize',10)


%% part B implicit method

% specify Reynolds number
Re=25; 

% set up the grids
Nx=51;Ny=21;   % no. of intervals
Lx=3;Ly=1;
dx=Lx/Nx; dy=Ly/Ny;
x=linspace(0,Lx,Nx+1);
y=linspace(0,Ly,Ny+1);
[x2,y2]=meshgrid(x,y);   

% specify the value of dt
dtim=0.02;  
t_f=0.2;
NMax=t_f/dtim;

% Assemble the A matrix to solve the Poisson equation Au=b
gamma = dx^2/dy^2;
A=Amatrix(Nx,Ny,dx^2,gamma);

% initial condition for the vorticity
vort = zeros(Nx-1,Ny-1);
tim=0;

% An array to record the kinetic energy at a certain point
energyarray=[];
uCimp =[]; tCimp =[];
figure

% time evoluion
for i=1:NMax

% solve possion eqn
stmfunc = solve_Poisson(vort,A,Nx,Ny);                                                             
[LHS,RHS] = assembleLHSRHS(Nx,Ny,stmfunc,vort,dx,dy,dt,Re,dtim);
LHS = sparse(LHS);RHS = sparse(RHS);
u0 = zeros((Nx-1)*(Ny-1),1);
resobj = 10^-4;
vort = SORmethod(LHS,RHS,u0,resobj);
vort = reshape(vort,Nx-1,Ny-1);

disp(['Finish the timestep ' num2str(i)])
tim = tim+dtim;

% update list
subplot(1,2,1)
[u,v]=get_uv(stmfunc,Nx,Ny,dx,dy);
uCimp = [uCimp; u(ic,jc)];
tCimp = [tCimp; tim];
% plot
quiver(x2,y2,u',v',5)
title(['At time step = ' num2str(i)])
set(gca,'FontSize',10)
xlabel('x');ylabel('y');ylim([0,1]);xlim([0,3]);
pause(0.2)
% drawnow
energyarray = [energyarray sqrt(u(end-1,end-1)^2+v(end-1,end-1)^2)];
subplot(1,2,2)
plot(dt:dt:dt*i,energyarray,'-*b')
title(['The kinetic energy at the grid point (' num2str(x(end-1)) ',' num2str(y(end-1)) ')'])
set(gca,'FontSize',10)
xlabel('time');ylabel('energy')

end

%% post-processing
plotImplicit(stmfunc,vort,Nx,Ny,dx,dy,x,y,tim);
% velocity plot
figure
subplot(1,2,1)
plot(tCimp,uCimp,'-b','LineWidth',1.5)    
title(['u velocity at point C(1.5,0.5) for Implicit Method'])    
set(gca,'FontSize',10)
xlabel('time')
ylabel('u signal')
% amplitude plot
Fs = 1/dtim;
L = length(uCimp);
Y = fft(uCimp);
P2 = abs(Y/L);
P1 = P2(1:(L+1)/2);
P1(2:end) = 2*P1(2:end);
ff = Fs*(0:(L-1)/2)/L;
subplot(1,2,2)
plot(ff,P1,"-x") 
title("Single-Sided Amplitude")
xlabel("f (Hz)")
ylabel("P1 (f)")
set(gca,'FontSize',10)

