function [ reflectedVector ] = reflect(vector, normal)
%UNTITLED2 Calculates the reflected vector using angle of incidence
%   equals angle of reflectance

reflectedVector = vector - normal * 2 * dot(vector, normal);

end