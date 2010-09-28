function [ R ] = fresnelCoefficient(vector, normal, n1, n2 )
%UNTITLED4 Computed fresnel equation for unpolarized light 
cosI = -dot(vector,normal);
sinISquared = 1 - cosI^2;
nSquared = (n1/n2)^2;
rootTerm = sqrt(1-nSquared*sinISquared);

rS = ((n1*cosI - n2*rootTerm) / (n1*cosI + n2*rootTerm))^2;
rP = ((n1*rootTerm - n2*cosI)/ (n1*rootTerm +n2*cosI))^2;

R = (rS + rP) / 2;
end