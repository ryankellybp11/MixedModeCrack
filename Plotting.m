% ======================================================================= %
% Solving the system of integral equations in the Konda/Erdogan paper     %
% Ryan Kelly                                                              %
% University of Houston Downtown, REU 2019                                %
% 7/12/2019                                                               %
%                                                                         %
% This code calculates the relative crack surface displacement and plots  %
% it on the domain for two combinations of delta and theta in comparison  %
% to a nonhomogeneous medium (delta = 0)                                  %
% ======================================================================= %
clear; clc; close all; format compact; format long g
set(0,'DefaultTextInterpreter','latex')


% Defining Constants
N = 40; % Number of collocation points (This must be an even number)
nu = 0.3; % Poisson's Ratio
An = zeros(N,2);
Bn = zeros(N,2);
kappa = 3 - 4*nu; % Plane strain
% kappa = (3 - nu)/(1 + nu); % Plane stress
UpBnd = 100;

for ij = 1:2
    if ij == 1
        delta = 2.5;
        theta = 0;
    else
        delta = 2.5;
        theta = pi/2;
    end
        
        B = delta.*cos(theta); % Beta
        G = delta.*sin(theta); % Gamma
        
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

        % Creating the coefficient matrices 
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
        An(:,ij) = X(1:N);
        Bn(:,ij) = X(N+1:end);
        
end

x = linspace(-.999,.999,200);
f1 = zeros(1,length(x));
f2 = zeros(1,length(x));
f3 = zeros(1,length(x));
n = 1:N;

for i = 1:length(x)
    f1(i) = sqrt(1 - x(i).^2);
    f2(i) = -2/(1+kappa)*sqrt(1 - x(i).^2).*sum(Bn(:,1)'.*u(x(i),0:N-1)./n); 
    f3(i) = -2/(1+kappa)*sqrt(1 - x(i).^2).*sum(Bn(:,2)'.*u(x(i),0:N-1)./n); 
end

plot(x,f1,'k.-')
hold on
plot(x,f2,'k-')
plot(x,f3,'k--')
xlabel('$x/a$')
ylabel('$\overline{v}$')
legend('\delta = 0','\theta = 0','\theta = \pi/2')

% for delta = 0.5
% axis([-1 1 0 1.2])
% yticks(0:0.4:1.2)

% for delta = 2.5
axis([-1 1 0 2.4])
yticks(0:0.4:2.4)