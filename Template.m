% ======================================================================= %
% Solving the system of integral equations in the Konda/Erdogan paper     %
% Ryan Kelly                                                              %
% University of Houston Downtown, REU 2019                                %
% 6/20/2019                                                               %
%                                                                         %
% This code was written as a template for T1 and T2, but it calculates all%
% 4 SIFs for a particular set of parameters. Be aware that small delta    %
% values can take a while to compute.                                     %
% ======================================================================= %
clear; clc; close all; format compact; format long g
set(0,'DefaultTextInterpreter','latex')

% Defining Constants
N = 30; % Number of collocation points (This must be an even number)
nu = 0.3; % Poisson's Ratio
delta = 0.25;
theta = 0*pi;
       
B = delta.*cos(theta); % Beta
G = delta.*sin(theta); % Gamma
kappa = 3 - 4*nu; % Plane strain
% kappa = (3 - nu)/(1 + nu); % Plane stress
UpBnd = 100;

% Creating the integration (s) and collocation (r) vectors
s = zeros(1,N+1); % s is a row vector
r = zeros(N,1); % r is a column vector
for ii = 1:N
    s(ii) = cos((2*ii-1)*pi/(N*2));
    r(ii) = cos(ii*pi/N);
end
s(N+1) = cos((2*(N+1)-1)*pi/(N*2));

% Loading Functions
p1 = -ones(N,1).*exp(-B.*r);
p2 = zeros(N,1); 

U = zeros(N,N); % Chebyshev Polynomials of the Second Kind
for jj = 1:N
    for ii = 1:N
        U(jj,ii) = u(r(jj),ii-1);
    end
end

T = zeros(N+1,N); % Chebyshev Polynomials of the First Kind
for m = 1:N+1
    for n = 1:N
        T(m,n) = t(s(m),n)/(N+1);
    end
end

L11 = zeros(N,N+1);
L12 = zeros(N,N+1);
L21 = zeros(N,N+1);
L22 = zeros(N,N+1);
for n = 1:N
    for m = 1:N+1
        fun1 = @(x) KK11(x,kappa,B,G,s(m),r(n));
        L11(n,m) = real(integral(fun1,0,UpBnd));

        fun2 = @(x) KK12(x,kappa,B,G,s(m),r(n));
        L12(n,m) = real(integral(fun2,0,UpBnd));

        fun3 = @(x) KK21(x,kappa,B,G,s(m),r(n));
        L21(n,m) = real(integral(fun3,0,UpBnd));

        fun4 = @(x) KK22(x,kappa,B,G,s(m),r(n));
        L22(n,m) = real(integral(fun4,0,UpBnd));
    end
end

% Creating the coefficient matrices from Eq. (37)
a = L11*T;
b = U + L12*T;
c = U + L21*T;
d = L22*T;

% Combining Systems
% I am putting both systems into a single matrix in order to solve it
ab = cat(2,a,b);
cd = cat(2,c,d);
M = cat(1,ab,cd);
p = cat(1,p1,p2);

X = M\p;

% Extracting An and Bn from result
An = X(1:N);
Bn = X(N+1:end);

% Calculating SIFs
k1a = -sum(Bn).*exp(B);
k2a = -sum(An).*exp(B);
k1na = 0;
k2na = 0;

for ii = 1:N
    if mod(ii,2) == 0
        k1na = k1na + Bn(ii);
        k2na = k2na + An(ii);
    else
        k1na = k1na - Bn(ii);
        k2na = k2na - An(ii);
    end
end
k1na = k1na.*exp(-B); % exp(-B) is because it is the left side of the crack
k2na = k2na.*exp(-B);

% Displaying results
fprintf('\n k1(a) = %.3f \n k1(-a) = %.3f \n k2(a) = %.3f \n k2(-a) = %.3f \n \n',k1a,k1na,k2a,k2na)
