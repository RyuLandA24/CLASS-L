function [qsat, P] = calcqsat(T,P)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

esat      = 0.611e3 * real(exp((17.2694 .* (T - 273.16)) ./ (T - 35.86)));
P         = P;
qsat      = 0.622 * esat ./ (P-(1-0.622).*esat);

end
% Unit is Pa, T is Kevin
