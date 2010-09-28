function [ refracted ] = refract( vector, normal, n1, n2)
%UNTITLED3 Computes refracted vector using Snell's law and Heckbert's method
cosI = -dot(vector,normal);
n = n1/n2;
refracted = n*vector + (n*cosI - sqrt(1 - n^2*(1-cosI^2)))*normal;

end