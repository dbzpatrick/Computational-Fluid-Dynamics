%% ME4233 Assignment 1 Part 2

% Dai Baizhou A0266515B
% Sam Marret A0252323R
% Aditiya Satish Nalini A0244516H


%% (a)

% C value
C = mean([26,025,024]);
Nx = 60; Ny =50;
Lx = 1; Ly = 1;
% initialize
I = Nx-1;
J = Ny-1;
IJ = I*J;
A = zeros(IJ);
dx = Lx/Nx; dy = Ly/Ny;
gamma = (4/3)*(dx/dy);
D = -2*(gamma^2+1);
% fill in D
i = IJ;
for j = IJ:-1:1
    A(j,i) = D;
    i = i-1;
end
% fill in gamma
i =IJ-Nx+1;
for j = IJ:-1:1
    A(j,i) = gamma^2;
    i = i-1;
    if i*j == 0
        break
    end
end
i =IJ;
for j = IJ-Nx+1:-1:1
    A(j,i) = gamma^2;
    i = i-1;
    if i*j == 0
        break
    end
end
% fill in 1
for i = 1:IJ-1
    if i/I == floor(i/I)
        A(i,i+1) =0;
        A(i+1,i) =0;
    else 
        A(i,i+1) =1;
        A(i+1,i) =1;
    end
end
% find b
% boundary values
B = zeros(Ny+1,Nx+1);
for i = 1:Nx+1
    x = (i-1)*dx;
    y = (i-1)*dy;
    B(end,i) = 0; %bottom
    B(1,i) = 0; %top
    B(i,1) = 0; %left
    B(i,end) = 0; %right
end
% g values and finalize b
idx = 0;
for j = Ny:-1:2
    for i =2:Nx     
        idx = idx+1;
        b(idx,1) = B(j-1,i)+B(j,i-1)+B(j,i)+B(j+1,i)+B(j,i+1);
        x = (i-1)*dx;
        y = (Ny-j+1)*dy;
        g = -C*(13/9)*pi^2*cos(2*pi*(2*x+y));
        b(idx,1) = 16*g*dx^2 -b(idx,1);
    end
end    
% Solve the linear problem Au=b 
A = sparse(A);
u_linear = A\b; 
N = size(A,1);
% linear solution plotting
u2D = reshape(u_linear,Nx-1,Ny-1);
u2D = [zeros(Nx-1,1) u2D zeros(Nx-1,1) ];
u2D= [zeros(1,Ny+1);u2D;2*ones(1,Ny+1) ];
figure
nexttile
surf(0:dx:dx*Nx,0:dy:dy*Ny,u2D')
xlabel('x');ylabel('y');zlabel('u')
title(['Solution to Au=b with linear method'])
set(gca,'FontSize',15)
xlim([0 1]);ylim([0 1])

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

%% QR decomposition solution
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
for omega = 1.4:0.1:1.9
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
% empty grid
x = dx:dx:dx*(Nx-1);
y = dy:dy:dy*(Ny-1);
[X, Y] = meshgrid(x, y);
% initial guess 1
u01 = ones(N,1);
u02 = C*cos(2*pi*(2*X+Y)); u02 = reshape(u02', [], 1);
u03 = 100*ones(N,1); 
% omega
omega = 1.5;
Asor = A;
figure
% SOR method with different initial guesses
for u0_sor = [u01,u02,u03]
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
    semilogy(0:k,resarray)
    set(gca,'FontSize',10)
    ylabel('Residual')
    xlabel('Iteration')
    title('Residual vs iteration number (different SOR initial guess)')
    hold on
end
legend('u01: vector of ones', 'u02: u = C cos(2pi(2x+y))','u03: constant vector of 100')

%% (e) change boundary conditions

% set up
Lx = 1; Ly = 1;
C = mean([26,025,024]);

% QR with different N
Nygrid = linspace(10,70,3);
RL2 = zeros(length(Nygrid),1);
h = zeros(length(Nygrid),1);
for i = 1:length(Nygrid)
    % initialize
    Ny = Nygrid(i);
    Nx = Ny+10;
    dx = Lx/Nx;
    dy = Ly/Ny;
    [Aqr,b] = findNewAb(Nx,Ny);

    %
    % find solution from QR    
    N = size(Aqr,1);    
    Q = zeros(N,N);
    % v1=u1
    Q(:,1)=Aqr(:,1);
    Q(:,1) = Q(:,1) /   sqrt( Q(:,1)'*Q(:,1) );
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
    %

    %Creating the analytical solution matrix
    idx = 1;
    u_analytical = zeros((Nx-1)*(Ny-1), 1); %(Nx2-1)*(Ny2-1) x 1 matrix
    for j = 1:Ny-1
        for i = 1:Nx-1
            idx = i + (j-1)*(Nx-1);
            x = dx*i;
            y = dy*j;
            u2 = C * cos(2* pi * (2 * x + y));
            u_analytical(idx) = u2;
        end
    end


    % find RL2
    RL2(i,1) = sqrt(sum((ut_qr-u_analytical).^2)/(Nx*Ny));

    % find h
    h(i,1) = sqrt(dx*dy);

    % plot
    % QR decomposition solution
    figure
    u2D_qr = reshape(ut_qr,Nx-1,Ny-1);
    u2D_qr = [zeros(Nx-1,1) u2D_qr zeros(Nx-1,1) ];
    u2D_qr = [zeros(1,Ny+1);u2D_qr;2*ones(1,Ny+1) ];
    nexttile
    surf(0:dx:dx*Nx,0:dy:dy*Ny,u2D_qr')
    xlabel('x');ylabel('y');zlabel('u')
    title(['Solution to Au=b with QR decompositions'])
    set(gca,'FontSize',15)
    xlim([0 1]);ylim([0 1])
end

figure
loglog(h,RL2,'-*y');    
hold on 
grid on
xlabel('h')
ylabel('RMS')
title('2nd-order accuracy verification (RMS vs h)')
loglog(h,100*h.^2,'-*r');
legend('Numerical results','f = c*h^2')
