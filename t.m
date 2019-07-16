%===============================================================
% Numerical Integral  Equation / for Vardoulakis Paper 
% Youn-Sha Chan
% 5/24/98      
%==============================================================

% This function computes the nth order Chebyshev polynomial of the first
% kind at x (either input may be a vector)

function out1 = T(x,n)
  out1 = cos(n*acos(x));