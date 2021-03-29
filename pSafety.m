%
% Script for computing p-safety /Rafael Wisniewski'21
% It uses Yalmit with sdpt3-solver
%
%
%%Initilization of the algorithm
clear 
x = sdpvar(2,1); %the states
%sdpvar p; %the probability to be found
%Unfortunately Yalmit is not mature to run full fledged optimisation
%therefore we have found p by multiple trials. Nonetheless, the code is
%written for "tru optimisation"; jest remore the comments of the command
%"sdpvar p" and comment the next line
p = 0.28;
D = 6; % the degree of all involved SOS polynomials
%%We use 8 unknown (sum of squares) polynomials
[s1,coefs1] = polynomial(x,D);
[s2,coefs2] = polynomial(x,D);
[s3,coefs3] = polynomial(x,D);
[s4,coefs4] = polynomial(x,D);
[s5,coefs5] = polynomial(x,D);
[s6,coefs6] = polynomial(x,D);
[s7,coefs7] = polynomial(x,D);
[s8,coefs8] = polynomial(x,D);

%%The switched diffusion system is considered with two subsystems and
%exponencial switched defined by  two transition rates
%The subsystems:
%   dx = f_1 dt + sigma_1 dW
%   dx = f_2 dt + sigma_2 dW
%The transition rates lambda_{12} and lambda_{21}.


sigma1 = .5*eye(2); % the diffusion matrix of System #1
sigma2 = .5*eye(2); % the diffusion matrix of System #2

f1 = [1; 1.4]; %drift of System #1
f2 = [1.4; 1]; %drift of System #1

lambda12 = 10;
lambda21 = 10;

%The difinition of three involved sets
%The state space  S = [S(i,x) \geq 0, i = 1,2] 
S1 = 10^2 - x(1)^2 - x(2)^2;
S2 = 10^2 - x(1)^2 - x(2)^2;
%The initial set A = [A(i,x) \geq 0, i = 1,2]
A1 = 1 - x(1)^2 - x(2)^2;
A2 = 1 - x(1)^2 - x(2)^2;
%The forbidden state  U = [U(i,x) \geq 0, i = 1,2]
U1 = 1 - (x(1)-5)^2 - (x(2)-5)^2;
U2 = 1 - (x(1)-5)^2 - (x(2)-5)^2;


%The two test functions for Subsystem#1 and Subsystem#2 h = (h_1,h_2)
[h1,coefh1] = polynomial(x,D);
[h2,coefh2] = polynomial(x,D);

%The infinitesmall generator consists of two components
Lh1 = 0.5 * trace(sigma1*sigma1'*hessian(h1,x)) + jacobian(h1,x)*f1 + lambda12*(h2-h1);
Lh2 = 0.5 * trace(sigma2*sigma2'*hessian(h2,x)) + jacobian(h2,x)*f2 + lambda21*(h1-h2);
    %Pure diffusion processes for the test
    %Lh1 = 0.5 * trace(sigma1*sigma1'*hessian(h1,x)) + jacobian(h1,x)*f1;
    %Lh2 = 0.5 * trace(sigma2*sigma2'*hessian(h2,x)) + jacobian(h2,x)*f2;
    %Pure dynamical (deterministic) systems for the test
    %Lh1 = jacobian(h1,x)*f1;
    %Lh2 = jacobian(h2,x)*f2;


%Constraints of the optimisation
%We use eight sos: s1, s2, s3, s4, s5, s6, s7, s8
constr = [sos(s1);sos(s2);sos(s3);sos(s4);sos(s5);sos(s6);sos(s7);sos(s8)];
%h >= 0 on the state space S 
constr = [constr; sos(h1 - s1*S1)];
constr = [constr; sos(h2 - s2*S2)];

%The infnitesmall generator -Lh >= 0 on the state space S 
constr = [constr; sos(-Lh1 - s3*S1)];
constr = [constr; sos(-Lh2 - s4*S2)];

% p - h >= 0 on A
constr = [constr; sos(p-h1 - s5*A1)];
constr = [constr; sos(p-h2 - s6*A2)];

% h - 1 >= 0 on U
constr = [constr; sos(h1-1 - s7*U1)];
constr = [constr; sos(h2-1 - s8*U2)];


%SOLUTION
params = [coefs1;coefs2;coefs3;coefs4;coefs5;coefs6;coefs7;coefs8;coefh1;coefh2];
%options = sdpsettings('solver','mosek');
options = sdpsettings('solver','sdpt3');
solvesos(constr,p,options,params);
      
 
 
 
