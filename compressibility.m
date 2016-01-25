% Hydrogen Compressibility Factor (Z) Equation
%
% NIST Revised Standardized Equation for Hydrogen Gas (EW Lemmon)
% Densitites for Fuel Consumption Apllications:  Volume 113, Number 6,
% November - December 2008
%
% Accuracy and Valid Input
%         0.01 % from 220 K to 1000 K with pressures up to 70 MPa
%  within 0.01 % from 255 K to1000 K with pressures to 120 MPa
%  within 0.1 % from 200 K to 1000 K up to 200 MPa.

% Inputs:
% temperature [ºC]
% pressure [MPa]
%
% Output:
% Z (unitless factor)

function [Z] = compressibility(temperature, pressure)

a = [0.05888460 -0.06136111 -0.002650473 0.002731125 0.001802374 -0.001150707 0.9588528*10^-4 -0.1109040*10^-6 0.1264403*10^-9];
b = [1.325 1.87 2.5 2.8 2.938 3.14 3.37 3.75 4.0];
c = [1.0 1.0 2.0 2.0 2.42 2.63 3.0 4.0 5.0];

sub = zeros(9,1);

for zz = 1:9;
    sub(zz) = a(zz)*(double((100)/(temperature+273.15)))^b(zz)*(pressure/1)^c(zz);
end

Z = 1 + sum(sub);

