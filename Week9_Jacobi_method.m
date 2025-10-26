clear

% initialisation
A=[8 2 0
  3 -5 7
  -2 1 9];

b=[12 14 27]';

u0=[3 2 1]';          % the initial guess

N=size(A,1);

D = diag(diag(A));    % get the diagnoal part of A in Jacobi method
R = A - D;            % get the remaining part of A

k=0;                  % initialise iteration number
resarray = [];        % initialise the array to record the residual
disp(R*u0)
disp(D\(b-R*u0))
%%
while 1
    
    u1 = D\(b-R*u0);                  % Jacobi method
    disp(u1)
    
    residual = norm(u1-u0);   % norm(b-A*u0)        % calculate the norm
    resarray = [resarray residual];   % record the residual
    
    if residual < 10^-12             % exiting condition
        break
    end
    
    u0 = u1;                          % update u0 with u1, so that in next iteration, u0 is u1
    
    k = k+1;                          % iteration number adds one
    disp(['We at the iteration ' num2str(k) ' with the residual= ' num2str(residual)])
end

figure
semilogy(0:k,resarray,'-*b')
set(gca,'FontSize',40)
ylabel('Residual')
xlabel('Iteration')


