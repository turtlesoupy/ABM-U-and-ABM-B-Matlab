function [reflectance, transmittance, absorptance] = ABM(azimuthalI, polarI, interfaceArray, nSamples)
    startState    = 1;
    endState      = length(interfaceArray) + 1;
    absorbedState = length(interfaceArray) + 2;
    
    numTransmitted = 0;
    numReflected   = 0;
    numAbsorbed    = 0;
    
    %vector version
    [x,y,z] = sph2cart(azimuthalI, polarI, 1);
    directions = zeros(3, nSamples);
    directions(1,:) = x;
    directions(2,:) = y;
    directions(3,:) = z;
    %Directions already normalized
    states = ones(1, nSamples);
    stateInterfaces = repmat(abm_interface(), 1, length(interfaceArray) + 2);
    stateInterfaces(1:length(interfaceArray)) = interfaceArray;
    
    %Prellocation
    numVectors = nSamples;
    normals = zeros(3, numVectors);
    interfaces(1:numVectors) = abm_interface();
    n1s = zeros(1, numVectors);
    n2s = zeros(1, numVectors);
    perturbanceReflects = zeros(1, numVectors);
    perturbanceRefracts = zeros(1, numVectors);
    refractStates = zeros(1, numVectors);
    freeLengths = zeros(1, numVectors);
    R = zeros(1, numVectors);
    
    vectorizedStep(1);
    loops = 0;
    while(~all( (states == endState) | (states == startState) | (states == absorbedState)))
        vectorizedStep(0);
        loops = loops + 1;
    end
    loops
    reflectance   = sum(states == startState)  ./ nSamples;
    transmittance = sum(states == endState) ./ nSamples;
    absorptance   = sum(states == absorbedState) ./ nSamples;
    
    function vectorizedStep(first)       
        interfaces = stateInterfaces(states);
        isBelow = directions(3,:) >= 0; %counterintuitive but true
        isAbove = ~isBelow;
        belowInterfaces = interfaces(isBelow);
        aboveInterfaces = interfaces(~isBelow);
        
        
        if ~first
            previouslyReflected = (states == 1);
        else
            previouslyReflected = (states == -1);
        end
        previouslyAbsorbed = (states == absorbedState);
        previouslyTransmitted = (states == endState);
        
        %Defaults - particle going from above to below
        normals(3, isAbove) = 1;
        n1s(isAbove)     = [aboveInterfaces.n1];
        n2s(isAbove)    = [aboveInterfaces.n2];
        perturbanceReflects(isAbove) = [aboveInterfaces.perturbanceDownTop];
        perturbanceRefracts(isAbove) = [aboveInterfaces.perturbanceDownBottom];
        refractStates(isAbove)       = states(isAbove) + 1;
        
        %Particle going from below to above
        n1s(isBelow)   = [belowInterfaces.n2];
        n2s(isBelow)   = [belowInterfaces.n1];
        normals(3,isBelow) = -1;
        perturbanceReflects(isBelow) = [belowInterfaces.perturbanceUpBottom];
        perturbanceRefracts(isBelow) = [belowInterfaces.perturbanceUpTop];
        
        refractStates(isBelow) = states(isBelow) - 1;
        
        freeLengths(:) = freePathLength(directions, normals, [interfaces.absorptionCoefficient]);
        isAbsorbed = freeLengths < [interfaces.thickness];
        
        R(:) = fresnelCoefficient(directions, normals, n1s, n2s);
        
        shouldRefract = rand(1,numVectors) >= R;
        shouldRefractMat = repmat(shouldRefract, 3, 1);
        
        states(shouldRefract)     = refractStates(shouldRefract);
           
        reflects      = reflect(directions, normals);
        directions(:,:) = refract(directions, normals, n1s, n2s);
        directions(~shouldRefractMat) = reflects(~shouldRefractMat);
        
        
        perturbances = perturbanceReflects;
        perturbances(shouldRefract) = perturbanceRefracts(shouldRefract);
        directions = brakkeScattering(directions, perturbances);
        states(isAbsorbed) = absorbedState;
        states(previouslyReflected) = 1;
        states(previouslyAbsorbed) = absorbedState;
        states(previouslyTransmitted) = endState;
    end

    function [pathLength] = freePathLength(vectors, normals, absorptionCoefficients)
        rands  = rand(1,length(absorptionCoefficients));
        cosI   = -dot(vectors,normals,1);
        pathLength = -(1./absorptionCoefficients).*log(rands).*cosI;
    end

    function [reflectedVector] = reflect(vector, normal)
        % Calculates the reflected vector using angle of incidence
        % equals angle of reflectance
        normalCoeff =  2 .* dot(vector, normal,1);
        reflectedVector = vector - normal .* repmat(normalCoeff,3,1);
    end

    function [refracted] = refract( vector, normal, n1, n2)
        %Computes refracted vector using Snell's law and Heckbert's method
        cosI = -dot(vector,normal,1);
        n = n1./n2;
        normalCoeff = (n.*cosI - sqrt(1 - n.^2.*(1-cosI.^2)));
        %Note: not accounting for total internal reflection on purpose
        %(handled in fresnel)
        refracted = repmat(n,3,1).*vector + repmat(normalCoeff,3,1).*normal;
    end

    function [R] = fresnelCoefficient(vector, normal, n1, n2 )
        %Computed fresnel equation for unpolarized light 
        cosI = -dot(vector,normal,1);
        sinISquared = 1 - cosI.^2;
        nSquared = (n1./n2).^2;
        rootTerm = 1-nSquared.*sinISquared;
        rootedTerm = sqrt(1-nSquared.*sinISquared);

        rS = ((n1.*cosI - n2.*rootedTerm) ./ (n1.*cosI + n2.*rootedTerm)).^2;
        rP = ((n1.*rootedTerm - n2.*cosI)./ (n1.*rootedTerm +n2.*cosI)).^2;
        
        R = (rS + rP) ./ 2;
        
        R(rootTerm < 0) = 1; %Total internal reflection
    end

    function [ perturbed ] = brakkeScattering(vector, delta)
        % Perturbs a vector according to exponentiated cosine distrubtion,
        % as described by Brakke et al. (1989)
        %assumes unit vectors
        theta = acos(vector(3, :));
        phi   = atan2(vector(2, :), vector(1,:));
        
        %[theta, phi, r] = cart2sph(vector(1), vector(2), vector(3));
        alpha = acos(rand(1,length(delta)).^(1./(delta + 1)));
        beta = 2.*pi.*rand(1,length(delta));
        perturbed = zeros(3,length(delta));
        
        newTheta = alpha + theta;
        newPhi   = phi + beta;
        st = sin(newTheta);
        sp = sin(newPhi);
        ct = cos(newTheta);
        cp = cos(newPhi);
        
        perturbed(1,:) = st .* cp;
        perturbed(2,:) = st .* sp;
        perturbed(3,:) = ct;
        
        %[perturbed(1), perturbed(2), perturbed(3)] = sph2cart(theta + beta, phi+alpha,r);
    end
end

