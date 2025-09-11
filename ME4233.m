%% assignment 1 ME4233

% matrix A
A = [-2 3/2 0 0; 1/2 -2 3/2 0; 0 1/2 -2 3/2; 0 0 1/2 -2];
[E, D] = eig(A);

disp('D');
disp(D);

disp('E');
disp(E);

%% assignment 2

% non-homogenous term f
f = [0.5;0;0;3/2];

% non-homogenous F
F = inv(E)*f;
disp('F');
disp(F);

% solve c1 to c4
syms c1 c2 c3 c4
eq1 = 1 == c1-(-1.4954/-3.4013);
eq2 = 1 == c2-(3.2966/-2.5352);
eq3 = 1 == c3-(-1.7004/-0.5987);
eq4 = 1 == c4-(2.8991/-1.4648); 
sol = solve([eq1, eq2, eq3, eq4], [c1, c2, c3, c4]); 
c1 = double(sol.c1);
c2 = double(sol.c2); 
c3 = double(sol.c3);
c4 = double(sol.c4);

% eigenvalues
ld1 = D(1,1); ld2 = D(2,2); ld3 = D(3,3); ld4 = D(4,4);

% matrix U
syms t
U = [c1*exp(ld1*t)-(F(1)/ld1);c2*exp(ld2*t)-(F(2)/ld2);c3*exp(ld3*t)-(F(3)/ld3);c4*exp(ld4*t)-(F(4)/ld4)];

% small u
u = E*U;
u = vpa(u,3);
disp('u');
disp(u);

% solve complimentery soln
syms sigma h ld
sol = solve(sigma^2-2*h*ld*sigma-1==0,sigma);
disp('complimentary solution:')
disp(sol)

% solve b1 b2
syms b1 b2
eq1 = 1== b1*(ld1+sqrt(ld1^2+1))+b2*(ld1-sqrt(ld1^2+1))+F(1)/ld1;
eq2 = 1== b1*(ld2+sqrt(ld2^2+1))+b2*(ld2-sqrt(ld2^2+1))+F(2)/ld2;
sol = solve([eq1,eq1], [b1,b2]);
b1 = vpa(sol.b1,4);
b2 = vpa(sol.b2,4);
disp('b1 and b2 solutions:');
disp([b1, b2]);
