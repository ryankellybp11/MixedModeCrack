%===============================================================
% Numerical Integral  Equation / for Vardoulakis Paper 
% Youn-Sha Chan
% 5/24/98      
%==============================================================


%This matlab program mainly defines the Chebyshev poly (in terms of xi)
%of 1st kind.


function out1 = T(x,n)
  out1 = cos(n*acos(x));


