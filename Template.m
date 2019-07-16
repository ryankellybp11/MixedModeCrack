% ======================================================================= %
% Solving the system of integral equations in the Konda/Erdogan paper     %
% Ryan Kelly                                                              %
% University of Houston Downtown, REU 2019                                %
% 6/20/2019                                                               %
% ======================================================================= %
clear; clc; close all; format compact; format long g
set(0,'DefaultTextInterpreter','latex')


% Defining Constants
N = 30; % Number of collocation points (This must be an even number)
nu = 0.3; % Poisson's Ratio
% FINAL = zeros(4,6);
% count = 0;
delta = 0.25;
theta = 0*pi;

% for delta = [0.1 0.25 0.5 1 2.5 5] % Delta, nonhomogeneity
%     count = count + 1;
    warning('off')
%     for theta = [0 pi.*0.5] % Angle of crack
        
        B = delta.*cos(theta); % Beta
        G = delta.*sin(theta); % Gamma
        kappa = 3 - 4*nu; % Plane strain
        % kappa = (3 - nu)/(1 + nu); % Plane stress
        LoBnd = -100;
        UpBnd = 100;

        % Creating the integration (s) and collocation (r) vectors
        s = zeros(1,N+1); % s is a row vector
        r = zeros(N,1); % r is a column vector
        for ii = 1:N
            s(ii) = cos((2*ii-1)*pi/(N*2));
            r(ii) = cos(ii*pi/N);
        end
        s(N+1) = cos((2*(N+1)-1)*pi/(N*2));
%         for n = 1:N
%             r(n) = cos(((n)*pi)/(N));
%         end
%         for n = 1:N+1
%             s(n) = cos((2*n+1)*pi/(2*(N+1)));
%         end



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
        % x = 0.001:0.001:100;
        for n = 1:N
            for m = 1:N+1
                fun1 = @(x) KK11(x,kappa,B,G,s(m),r(n));
                L11(n,m) = real(integral(fun1,0,UpBnd));

                fun2 = @(x) KK12(x,kappa,B,G,s(m),r(n));
                L12(n,m) = real(integral(fun2,0,UpBnd));

                fun3 = @(x) KK21(x,kappa,B,G,s(m),r(n));
                L21(n,m) = real(integral(fun3,0,UpBnd)); % for some reason, Matlab won't integrate this function from 0 for particular values of r and s

                fun4 = @(x) KK22(x,kappa,B,G,s(m),r(n));
                L22(n,m) = real(integral(fun4,0,UpBnd));
            end
        end

        % Creating the coefficient matrices from Eq. (37)

        a = L11*T;
        b = U + L12*T;
        c = U + L21*T;
        d = L22*T;
        % for n = 1:N
        %     for m = 1:N+1
        %         a(n,m) = 1/N.*(L11(n,m).*t(s(m),m));
        %         b(n,m) = u(r(n),m-1) + 1/N.*(L12(n,m).*t(s(m),m));
        %         c(n,m) = u(r(n),m-1) + 1/N.*(L21(n,m).*t(s(m),m));
        %         d(n,m) = 1/N.*(L22(n,m).*t(s(m),m));
        %     end
        % end

        % % Single-Valued Condition
        % Asv = zeros(1,N);
        % for n = 1:N
        %     
        %     f = @(s) t(s,n)./sqrt(1-s.^2);
        %     Asv(n) = integral(f,-1,1);   
        %     
        % end
        % Bsv = Asv;

        % Combining Systems
        % I am putting both systems into a single matrix in order to solve it
        ab = cat(2,a,b);
        cd = cat(2,c,d);
        % A_sv = cat(2,Asv,zeros(1,N));
        % B_sv = cat(2,zeros(1,N),Bsv);
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
        
%         if theta == 0 && delta ~= 5
%             FINAL(3,count) = k1a;
%             FINAL(4,count) = k1na;
%         elseif theta ~= 0
%             FINAL(1,count) = k1a;
%             FINAL(2,count) = k2a;
%         end
        
        % Displaying results
        fprintf('\n k1(a) = %.3f \n k1(-a) = %.3f \n k2(a) = %.3f \n k2(-a) = %.3f \n \n',k1a,k1na,k2a,k2na)
%     end
% end


%{
%% Plotting Checks
% This is examining the difference between my kernels (KKjj) and the
% paper-based kernels (Kjj)

x = linspace(-1000,1000,2000);
f1 = @(x) K11(x,kappa,B,G,s(end),r(end));
f2 = @(x) K12(x,kappa,B,G,s(end),r(end));
f3 = @(x) K21(x,kappa,B,G,s(end),r(end));
f4 = @(x) K22(x,kappa,B,G,s(end),r(end));


subplot(2,2,1)
hold on
plot(x,real(f1(x)),'b')
plot(x,real(fun1(x)),'r-.')
legend('K11','KK11')
grid on

subplot(2,2,2)
hold on
plot(x,real(f2(x)),'b')
plot(x,real(fun2(x)),'r-.')
legend('KK12','K12')
grid on

subplot(2,2,3)
hold on
plot(x,real(f3(x)),'b')
plot(x,real(fun3(x)),'r-.')
legend('K21','KK21')
grid on

subplot(2,2,4)
hold on
plot(x,real(f4(x)),'b')
plot(x,real(fun4(x)),'r-.')
legend('K22','KK22')
grid on

%}