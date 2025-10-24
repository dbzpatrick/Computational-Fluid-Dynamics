clear

% initialisation
% using SOR to solve Au=b for this matrix may diverge
% A=[  8    2   0
%      3   -5   7
%      -2   1   9]; % Non-SPD

A=[  8   3   0
     3   5   2
     0   2   9];   % This A is SPD

b=[12
   14
   27];

N=size(A,1);       % get the size of A

u0=[0 2 1]';       % the initial guess

omega=1.5;    % 0<omega<2     % the omega parameter in SOR

% extract the lower triangular entries from A to form a lower triangular matrix
L = tril(A,-1);
D=diag(diag(A));   % get the diagonal part of A
U=A-L-D;           % get the remaining part of A

k=0;               % k is the iteration number, starting from 0
resarray=[];       % record the residual

while 1
    u1=(D+omega*L)\(omega*b-(omega*U+(omega-1)*D)*u0);  % SOR method
    
    res=norm(u0-u1);            % calculate the norm
    resarray=[resarray res];    % record the residual
    
    if res<10^-12               % exiting condition
        break
    end
    
    u0=u1;                      % update u0 with u1, so that in next iteration, u0 is u1
    
    k=k+1;                      % iteration number adds one
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])
end

figure
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',40)

% The matrix G and its eigenvalues
G=-(D+omega*L)\(omega*U+(omega-1)*D);
Lambda=eig(G);

%plot a circle
% theta=0:0.02:2*pi;
% x = 1*cos(theta); y=1*sin(theta);
% 
% figure
% plot([-1.5 1.5],[0 0],'-k')
% hold on
% plot([0 0],[-1.5 1.5],'-k')
% plot(x,y,'-b')
% axis equal
% xlabel('Re(\sigma)');ylabel('Im(\sigma)')
% set(gca,'FontSize',30)
