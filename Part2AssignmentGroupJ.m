%% ME4233 Assignment 1 Part 2

% Dai Baizhou A0266515B
% Sam Marret A0252323R
% Aditiya Satish Nalini A0244516H


%% (a)
% C value
C = mean([26,025,024]);
% initialize
Nx = 60; Ny = 50;
Lx = 1; Ly = 1;
dx = Lx/Nx; dy = Ly/Ny;
temp = 16*dx^2;
Gamma = temp/(9*dy^2); D = -2*(Gamma+1);
% empty grid
x = dx:dx:dx*(Nx-1);
y = dy:dy:dy*(Ny-1);
% form matrix A
aux1 = toeplitz([D 1 zeros(1,Nx-3)],[D 1 zeros(1,Nx-3)]');
aux2 = Gamma^2*eye(Nx-1);
A=[aux1 aux2 zeros(Nx-1,(Ny-3)*(Nx-1))];
for i=1:Ny-3
    A = [A
         zeros(Nx-1,(i-1)*(Nx-1)) aux2 aux1 aux2 zeros(Nx-1,(Ny-i-3)*(Nx-1))];
end
A=[A
   zeros(Nx-1,(Ny-3)*(Nx-1))  aux2 aux1];
A = sparse(A);
% Form the vector b in Au=b
[X, Y] = meshgrid(x, y);
g = -C*(13/9)*(pi^2)*cos(2*pi*(2*X + Y));
g_vec = reshape(g', [], 1);
b = temp * g_vec;
% Boundary conditions for b 
bc=reshape([zeros(Nx-2,Ny-1);2*ones(1,Ny-1)],[],1);
b=b-bc;
% Solve the linear problem Au=b 
u_linear = A\b; 
N = size(A,1);

%% (b) 
% LU decomposition to solve for u
tic
Alu = A;
Lp = eye(N,N);  
% the loop to go from the 1st column to the second last column (N-1)
for i=1:N-1
    auxL=eye(N,N);    % auxiliary L matrix starts with an identity matrix
    
    auxL(i+1:N,i) = - Alu(i+1:N,i)/Alu(i,i);  % the equation to calculate the coefficients in the L matrix
    Alu(i+1:N,:) = Alu(i+1:N,:) + auxL(i+1:N,i)*Alu(i,:);  % updating the A matrix using the row operations. 
    Lp=auxL*Lp;   % save the auxiliary L matrix. 
end
U = Alu;
% get the intermediate vector yy=Lp*b;
y_lu = Lp*b;
% solve for u. Note that here I just use inv(U) directly. 
ut_lu=U\y_lu;
disp(ut_lu)
toc

%% QR decomposition
Aqr = A;
Q = zeros(N,N);
% v1=u1
Q(:,1)=Aqr(:,1);
Q(:,1) = Q(:,1) /   sqrt( Q(:,1)'*Q(:,1) ); % normalisation
% get v2,v3...,vj,vN in the orthogonal matrix
for j=2:N
    Proj=zeros(N,1);
    % to get vj, you first do projection of uj on v1--v_k---v_{j-1}
    for k=1:j-1
        Proj = Proj + ( Q(:,k)'*Aqr(:,j)  )/( Q(:,k)'*Q(:,k)   )*Q(:,k);
    end
    % and then subtract the projection from uj
    Q(:,j)  =  Aqr(:,j)  - Proj;
    % normalisation
    Q(:,j) = Q(:,j) /   sqrt( Q(:,j)'*Q(:,j) );
end
R=zeros(N,N);
% get R
for j=1:N
    % you do the projection of uj on e1,e2,...,e_j
    for k=1:j
        R(k,j) = Q(:,k)'*Aqr(:,j);
    end
end
% solve
y_qr=Q'*b;
ut_qr=R\y_qr;

%% plotting
% LU decomposition solution
u2D_lu = reshape(ut_lu,Nx-1,Ny-1);
u2D_lu = [zeros(Nx-1,1) u2D_lu zeros(Nx-1,1) ];
u2D_lu = [zeros(1,Ny+1);u2D_lu;2*ones(1,Ny+1) ];

figure
nexttile
surf(0:dx:dx*Nx,0:dy:dy*Ny,u2D_lu')
xlabel('x');ylabel('y');zlabel('u')
title(['Solution to Au=b with LU decompositions'])
set(gca,'FontSize',15)
xlim([0 1]);ylim([0 1])

% QR decomposition solution
u2D_qr = reshape(ut_qr,Nx-1,Ny-1);
u2D_qr = [zeros(Nx-1,1) u2D_qr zeros(Nx-1,1) ];
u2D_qr = [zeros(1,Ny+1);u2D_qr;2*ones(1,Ny+1) ];
nexttile
surf(0:dx:dx*Nx,0:dy:dy*Ny,u2D_qr')
xlabel('x');ylabel('y');zlabel('u')
title(['Solution to Au=b with QR decompositions'])
set(gca,'FontSize',15)
xlim([0 1]);ylim([0 1])

%% (c)
random = rand(N,1);
% Jacobi method
Ajac = A;
u0_jac = random;   % random kick start vector
D = diag(diag(Ajac));    % get the diagnoal part of A in Jacobi method
R = Ajac - D;            % get the remaining part of A
k=0;                  % initialise iteration number
resarray = [];        % initialise the array to record the residual
while 1    
    u1 = D\(b-R*u0_jac);                  % Jacobi method    
    residual = norm(u1-u0_jac);   % norm(b-A*u0)        % calculate the norm
    resarray = [resarray residual];   % record the residual    
    if residual < 10^-7            % exiting condition
        break
    end    
    u0_jac = u1;                          
    k = k+1;                          % iteration number adds one
end
figure
nexttile
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',10)
ylabel('Residual')
xlabel('Iteration')
title('Residual vs iteration number (Jacobi method)')

%% Gauss-Seidel method
u0_gs = random;
Ags = A;
% extract the lower triangular entries from A to form a lower triangular matrix
L = tril(Ags);
U=Ags-L;             % get the remaining part of A
k=0;               % k is the iteration number, starting from 0
resarray=[];       % record the residual
while 1    
    u1=L\(b-U*u0_gs);                 % Gauss-Seidel method     
    res=norm(u0_gs-u1);               % calculate the norm
    resarray=[resarray res];       % record the residual
    if res<10^-7                   % exiting condition
        break
    end    
    u0_gs=u1;                         % update u0 with u1, so that in next iteration, u0 is u1    
    k=k+1;                         % iteration number adds one
end
nexttile
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',10)
ylabel('Residual')
xlabel('Iteration')
title('Residual vs iteration number (GS method)')

%% SOR method
Asor = A;
u0_sor = random;
omega = 1.5;
% extract the lower triangular entries from A to form a lower triangular matrix
L = tril(Asor,-1);
D=diag(diag(Asor));   % get the diagonal part of A
U=Asor-L-D;           % get the remaining part of A
k=0;               % k is the iteration number, starting from 0
resarray=[];       % record the residual
while 1
    u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0_sor);  % SOR method    
    res=norm(u0_sor-u1);            % calculate the norm
    resarray=[resarray res];    % record the residual    
    if res<10^-7              % exiting condition
        break
    end    
    u0_sor=u1;                      % update u0 with u1    
    k=k+1;                      % iteration number adds one
end
nexttile
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',10)
ylabel('Residual')
xlabel('Iteration')
title('Residual vs iteration number (SOR method)')

%% different omega value for SOR
figure
for omega = 0.8:0.2:1.8
    Asor = A;
    u0_sor = random;
    L = tril(Asor,-1);
    D=diag(diag(Asor));   % get the diagonal part of A
    U=Asor-L-D;           % get the remaining part of A
    k=0;               % k is the iteration number, starting from 0
    resarray=[];       % record the residual
    while 1
        u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0_sor);  % SOR method
        res=norm(u0_sor-u1);            % calculate the norm
        resarray=[resarray res];    % record the residual
        if res<10^-7              % exiting condition
            break
        end
        u0_sor=u1;                      % update u0 with u1
        k=k+1;                      % iteration number adds one
    end
    semilogy(0:k,resarray,'DisplayName',['omega = ' num2str(omega)]);
    set(gca,'FontSize',10)
    ylabel('Residual')
    xlabel('Iteration')
    sgtitle('Residual vs iteration number (SOR method)')
    legend
    hold on
end

%% (d)
% different initial guess for SOR
% initial guess 1
u01 = rand(N,1);
u02 = C*cos(2*pi*(2*X+Y)); u02 = reshape(u02', [], 1);
u03 = zeros(N,1); u03(:,1) = 1; u03 = reshape(u03', [], 1);
% omega
omega = 1.5;
Asor = A;
% SOR method with different initial guesses
for u0 = {u01, u02, u03}
    u0_sor = u0;

    % extract the lower triangular entries from A to form a lower triangular matrix
    L = tril(Asor,-1);
    D=diag(diag(Asor));   % get the diagonal part of A
    U=Asor-L-D;           % get the remaining part of A
    k=0;               % k is the iteration number, starting from 0
    resarray=[];       % record the residual
    while 1
        u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0_sor);  % SOR method
        res=norm(u0_sor-u1);            % calculate the norm
        resarray=[resarray res];    % record the residual
        if res<10^-7              % exiting condition
            break
        end
        u0_sor=u1;                      % update u0 with u1
        k=k+1;                      % iteration number adds one
    end
    semilogy(0:k,resarray,'-*b')
    set(gca,'FontSize',10)
    ylabel('Residual')
    xlabel('Iteration')
    title('Residual vs iteration number (SOR method)')
    hold on
end