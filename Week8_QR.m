% Using the modified GS for QR decomposition to solve Au=b
clear
close all

% Initialisation
A=[1 3 2 
   0 -1.6 -1.6
   0 -1.2 3.8];

N=size(A,1);

Q = zeros(N,N);

% v1=u1
Q(:,1)=A(:,1);
Q(:,1) = Q(:,1) /   sqrt( Q(:,1)'*Q(:,1) ); % normalisation

% get v2,v3...,vj,vN in the orthogonal matrix
for j=2:N
    Proj=zeros(N,1);
    % to get vj, you first do projection of uj on v1--v_k---v_{j-1}
    for k=1:j-1
        Proj = Proj + ( Q(:,k)'*A(:,j)  )/( Q(:,k)'*Q(:,k)   )*Q(:,k);
    end
    % and then subtract the projection from uj
    Q(:,j)  =  A(:,j)  - Proj;
    % normalisation
    Q(:,j) = Q(:,j) /   sqrt( Q(:,j)'*Q(:,j) );

    disp()
end

R=zeros(N,N);

% get R
for j=1:N
    % you do the projection of uj on e1,e2,...,e_j
    for k=1:j
        R(k,j) = Q(:,k)'*A(:,j);
    end
end

b=[13 -8 9]';

y=Q'*b;
u=R\y;


