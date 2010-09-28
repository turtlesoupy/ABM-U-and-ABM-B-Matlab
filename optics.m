function [ refracted ] = refract( vector, normal, n  )
%UNTITLED3 Computes refracted vector using Snell's law and Heckbert's method
cosI = -dot(vector,normal);
refracted = n*vector + (n*cosI - sqrt(1 - n^2*(1-cosI^2)))*normal;

end

function [ reflectedVector ] = reflect(vector, normal )
%UNTITLED2 Calculates the reflected vector using angle of incidence
%   equals angle of reflectance

reflectedVector = vector - normal * 2 * dot(vector, normal);

end

function [ R ] = fresnelReflectanceUnpolarized(vector, normal, n1, n2 )
%UNTITLED4 Computed fresnel equation for unpolarized light 
cosI = -dot(vector,normal);
sinISquared = 1 - cosI^2;
nSquared = (n1/n2)^2;
rootTerm = sqrt(1-nSquared*sinISquared);

rS = ((n1*cosI - n2*rootTerm) / (n1*cosI + n2*rootTerm))^2;
rP = ((n1*rootTerm - n2*cosI)/ (n1*rootTerm +n2*cosI))^2;

R = (rS + rP) / 2;
end

function [ perturbed ] = brakkeScattering(vector, delta)
%UNTITLED5 Perturbs a vector according to exponentiated cosine distrubtion,
% as described by Brakke et al. (1989)
[theta, phi, r] = cart2sph(vector(1), vector(2), vector(3));
alpha = arccos(rand()^(1/(delta + 1)));
beta = 2*pi*rand();
pertubed = sph2cart(theta + alpha, phi+beta,r);

end