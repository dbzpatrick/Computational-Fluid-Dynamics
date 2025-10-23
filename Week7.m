% ME4233 Week 7 code for solving a Poisson equation
% You don't have to understand this code. Just execute it for fun. You can 
% for example change the values of M, N or bctype below to see the effect.
% We will explain how to solve a Poisson equation in the following weeks 
% step by step.

% clear
% close all

% Initialisation
% you can change the values of M and N to see the effect of grid number
M=40;  Lx=1;    dx=Lx/M; % Delta x
N=50;  Ly=0.5;  dy=Ly/N;
Gamma = dx/dy;   D = -2*(Gamma^2+1); % some auxiliary constants

% Grid or mesh
x = dx:dx:dx*(M-1);
y = dy:dy:dy*(N-1);

% Assemble the A matrix in Au=b
% the following lines (24-33) are difficult to understand. You don't have
% to understand them. You will learn a more intuitive way to assemble A
% later.
aux1 = toeplitz([D 1 zeros(1,M-3)],[D 1 zeros(1,M-3)]');
aux2 = Gamma^2*eye(M-1);

A=[aux1 aux2 zeros(M-1,(N-3)*(M-1))];
for i=1:N-3
    A = [A
         zeros(M-1,(i-1)*(M-1)) aux2 aux1 aux2 zeros(M-1,(N-i-3)*(M-1))];
end
A=[A
   zeros(M-1,(N-3)*(M-1))  aux2 aux1];

% Form the vector b in Au=b
b = 1*dx^2*ones((M-1)*(N-1),1);

% Boundary conditions included in b
% You can change the type of bc below. Type 1 is the one in the slides 
% (note that it suffers from discontinuity at the two corners). Type 2 is
% smoother. Note the difference of implementation for the two bc's.
bctype=2;
if bctype==1
% uL=0, uR=2, uT=0, uB=0
bc=reshape([zeros(M-2,N-1);2*ones(1,N-1)],[],1);
elseif bctype==2
% uL=0, uR=sin(y/Ly*pi), uT=0, uB=0
bc=reshape([zeros(M-2,N-1);sin(y/Ly*pi)],[],1);
end
b=b-bc;

% Solve the linear problem Au=b using matlab command backslash \, but we 
% will learn other methods to solve the linear algebraic equation in the 
% next weeks
u = A\b; % to solve Au=b

% Reshape u into 2D matrix for plotting
u2D = reshape(u,M-1,N-1);
u2D = [zeros(M-1,1) u2D zeros(M-1,1) ];
if bctype==1
u2D = [zeros(1,N+1);u2D;2*ones(1,N+1) ];
elseif bctype==2
u2D = [zeros(1,N+1);u2D;[0 sin(y/Ly*pi) 0] ];
end

set(0,'DefaultFigureWindowStyle','docked') 
figure(3)
surf(0:dx:dx*M,0:dy:dy*N,u2D')
xlabel('x');ylabel('y');zlabel('u')
title(['Solution to a Poisson equation with M=' num2str(M) ', N=' num2str(N) ', bctype=' num2str(bctype)])
set(gca,'FontSize',30)
xlim([0 1]);ylim([0 0.5])
