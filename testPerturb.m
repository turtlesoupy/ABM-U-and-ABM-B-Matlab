function [] = testPerturb()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



vector = [1 -2 -1];
vector = vector/ norm(vector)

xyLength = sqrt(vector(1)^2 + vector(2)^2);

if xyLength == 0
    if vector(1)  > 0 
        zAngle = pi/2;
    else 
        zAngle = -pi/2;
    end
else
    zAngle = acos(vector(2)/xyLength);
end

xAngle = acos(xyLength);

if vector(3) < 0
    xAngle = -xAngle;
end

if vector(1) > 0
    zAngle = -zAngle
end


returnX = [1 0 0 ; 0 cos(xAngle) -sin(xAngle); 0 sin(xAngle) cos(xAngle)];
returnZ = [cos(zAngle) -sin(zAngle) 0; sin(zAngle) cos(zAngle) 0; 0 0 1];
%returnX = [1 0 0; 0 cos(azimuthal) sin(azimuthal); 0 -sin(azimuthal) cos(azimuthal)];

perturbVector = [0 1 0]';


perturbVector

returnZ * returnX *  perturbVector

end

