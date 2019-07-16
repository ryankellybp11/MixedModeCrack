%===============================================================
% Numerical Integral  Equation / for Vardoulakis Paper 
% Youn-Sha Chan
% 5/24/98      
%==============================================================


%This matlab program mainly defines the Chebyshev poly (in terms of xi)
%of 2nd kind.


function out1 = u(x,n)
  out1 = sin((n+1).*acos(x))./sin(acos(x));
end

