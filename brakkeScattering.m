function [ perturbed ] = brakkeScattering(vector, delta)
%UNTITLED5 Perturbs a vector according to exponentiated cosine distrubtion,
% as described by Brakke et al. (1989)
[theta, phi, r] = cart2sph(vector(1), vector(2), vector(3));
alpha = acos(rand()^(1/(delta + 1)));
beta = 2*pi*rand();
perturbed = zeros(1,3);
[perturbed(1), perturbed(2), perturbed(3)] = sph2cart(theta + beta, phi+alpha,r);
end