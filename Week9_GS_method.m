clear

% initialisation
A=[  8    2   0
     3   -5   7
     -2   1   9];
 
b=[12
   14
   27];

u0=[3 2 1]';       % the initial guess

N=size(A,1);       % get the size of A

% extract the lower triangular entries from A to form a lower triangular matrix
L = tril(A);
U=A-L;             % get the remaining part of A

k=0;               % k is the iteration number, starting from 0
resarray=[];       % record the residual

while 1
    
    u1=L\(b-U*u0);                 % Gauss-Seidel method
     
    res=norm(u0-u1);               % calculate the norm
    resarray=[resarray res];       % record the residual
    
    if res<10^-12                  % exiting condition
        break
    end
    
    u0=u1;                         % update u0 with u1, so that in next iteration, u0 is u1
    
    k=k+1;                         % iteration number adds one
    disp(['Finish iteration ' num2str(k) '. res = ' num2str(res)])
end

figure
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',40)