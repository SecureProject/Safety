%
%
% Script for computing weak p-safety /Rafael Wisniewski'21
% 
% The script uses bensolve tools
% Before running the script go to http://bensolve.org/tools/download.html.
%
clear all
% Degree of Bernstein basis
d = 30;
% Initilize a) the moments of the initial measure and the Bernstein
% b) coeficients of the indicator function I_[0.8, 1]
% c) The matrix with coefficients with he derivative
Y0 = zeros(d+1,1);
bq = zeros(d+1,1);
DF = zeros(d,d+1);
%Indicator function of I_[a1, b1]
a1 = 0.2; b1 = 1;
Ind = @(t) (t>=a1) & (t<=b1);
%The initial state is uniformly distributed on the interval [a2, b2]
a2 = 0; b2 = 0.1;
% m-th Bernstein basis
for m = 0:d
    Binom = factorial(d)/factorial(m)/factorial(d-m);
    Bern_m = @(x) Binom.*x.^m.*(1 - x).^(d-m);
    %%Computation of the Bernetein moments of uniform distribution on the interval [0, 0.2] 
    Y0(m+1) = 1/(b2-a2)*integral(Bern_m,a2,b2);
    % Compute Berstein approximation of the indicator function I_[0.8, 1] 
    bq(m+1) = 1.0*Ind(m/d);
    DF(m+1,m+1) = -1;
    DF(m+1,m+2)= 1;
end;

%Infinitesmall genertor of one-dimensional Brownina motion
D1F = d*DF(1:d,1:d+1); 
D2F = (d-1)*DF(1:d-1,1:d); 
DDF = 0.5*D2F*D1F;

% H-representation (hyperplane) in BenSolve
% P = [a  <= BX <= b, l <= x <= u]
clear rep;
rep.B= [-DDF;
        eye(d+1)];    
    
rep.a=zeros(length(rep.B),1);
P=polyh(rep,'h');
P=eval(P);
Rh = hrep(P);

% V-representation (vertices)
%P ={x ∈Rn|x =Vλ+Dμ+Lη, λ≥0, μ≥0, e⊤λ=1}
Rv = vrep(P);
% For a cone V and L are zero
% The dual cone to [BX >= 0] is [DX >= 0] 


% Linear optimisation for computing weak p-safety for one-dimensional Wiener process
%All the enequalities are of the form AX \leq b
% The definition of the enequality constraints
% B defines the dual cone [BX \geq 0]
B = Rv.D';


%combines constrainst AY <= b
A = [B;
    -eye(d+1);
    eye(d+1)];

b = [B*Y0; zeros(d+1,1); ones(d+1,1)];



% The definition of equality constrains
Aeq = ones(1,d+1);
beq = 1;
% The definition of the (linear) objective function
f = bq;

% The solution of minimisation subject to the constraints Ax <= b and Aeq x = beq
% 
options = mskoptimset('');
options = mskoptimset(options,'Diagnostics','on');
[x,fval] = linprog(-f,A,b,Aeq,beq,[],[],options)
-fval