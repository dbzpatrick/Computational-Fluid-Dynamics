% Using LU decomposition to solve Au=b
clear
close all

% Initialisation
% A0=[3 2 1
%     6 2 3
%    -3 -3 8.5];
% 
% N=size(A0,1);
N=10;
A0 = rand(N,N);

% Lp matrix starts with an identity matrix. Initialize the U matrix (which
% is redundant).
Lp=eye(N,N);

% make a copy of A in U and we do the row operation on U
A=A0;

% the loop to go from the 1st column to the second last column (N-1)
for i=1:N-1
    auxL=eye(N,N);    % auxiliary L matrix starts with an identity matrix
    
    auxL(i+1:N,i) = - A(i+1:N,i)/A(i,i);  % the equation to calculate the coefficients in the L matrix
    A(i+1:N,:) = A(i+1:N,:) + auxL(i+1:N,i)*A(i,:);  % updating the A matrix using the row operations. 
                                                     % You are using U(j,:) row to do these row operations.
    
    Lp=auxL*Lp;   % save the auxiliary L matrix. 
                  % Remember that the last operation will left-multiply the previous operations. 
                  % So auxL appears on the left in auxL*Lp.
end
U = A;
% We got A=LU. Now we can use it to solve Au=b.

% b vector. When boundary conditions change or the source terms change, you
% can still use the LU decomposition of A above to solve Au=b (for
% different b).
b=[11 18 5]';

% get the intermediate vector y
y=Lp*b;

% solve for u. Note that here I just use inv(U) directly. 
u=U\y;









% or you can implement the back elimination as
for i=N:-1:1
    R = 0;
    for j = i+1:N
        R = R + U(i,j)*u1(j);
    end
    u1(i) = (y(i) - R )/U(i,i);
end









