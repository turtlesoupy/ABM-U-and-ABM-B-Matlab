function [reflectance, transmittance, absorptance] = ABM(azimuthalI, polarI, interfaceArray, nSamples)
    stepFunction  = @step;    
    startState    = 1;
    endState      = length(interfaceArray) + 1;
    absorbedState = length(interfaceArray) + 2;
    endStates = zeros(1, nSamples);
    [x,y,z] = sph2cart(azimuthalI, polarI, 1);
    startDirection = [x,y,z];
    
    startSplits = find([interfaceArray.splitThicknessIndex]);
    endSplits   = [interfaceArray(startSplits).splitThicknessIndex];
    numSplits = length(startSplits);
    baseThicknesses = [interfaceArray.thickness];
    parfor i = 1:nSamples
        direction = startDirection;
        state = startState;

        %ABM-U splits some of the thickness
        thicknesses = baseThicknesses;
        if numSplits > 0
            rands = rand(1,numSplits);
            thicknesses(startSplits) = thicknesses(startSplits) .* rands; 
            thicknesses(endSplits)   = thicknesses(endSplits)   .* (1-rands);
        end
        
        [state, direction] = stepFunction(state, direction, thicknesses);
        while(state ~= endState && state ~= 1 && state ~= absorbedState)
            [state, direction] = stepFunction(state, direction, thicknesses);
        end
        
        endStates(i) = state;
        
    end
    
    reflectance   = sum(endStates == startState) / nSamples;
    transmittance = sum(endStates == endState) / nSamples;
    absorptance   = sum(endStates == absorbedState) / nSamples;
    
    function [outState, outVector] = step(state, vector, interfaceThicknesses)
        interface = interfaceArray(state);
        thickness = interfaceThicknesses(state);
        if vector(3) < 0 
            normal = [0,0,1];
            n1 = interface.n1;
            n2 = interface.n2;
            perturbanceReflect = interface.perturbanceDownTop;
            perturbanceRefract = interface.perturbanceDownBottom;
            refractState = state + 1;
        else
            normal = [0,0,-1];
            n1 = interface.n2;
            n2 = interface.n1;
            perturbanceReflect = interface.perturbanceUpBottom;
            perturbanceRefract = interface.perturbanceUpTop;
            refractState = state - 1;
        end

        normalAngle = -dot(vector, normal);

        %Check that we aren't absorbed here        
        if thickness > 0 && ...
                freePathLength(vector, normal, normalAngle, interface.absorptionCoefficient) < thickness
            outState = absorbedState;
            outVector = [0 0 0];
        else
            R = fresnelCoefficient(vector, normal, normalAngle, n1, n2);
            if rand() < R
                outState = state;
                outVector = reflect(vector, normal, normalAngle);
                if ~isinf(perturbanceRefract)
                    outVector = brakkeScattering(outVector, perturbanceReflect);
                end
            else
                outState  = refractState;
                outVector = refract(vector, normal, normalAngle, n1, n2);
                if ~isinf(perturbanceRefract)
                    outVector = brakkeScattering(outVector, perturbanceRefract);
                end
            end
            
        end                 
    end

    function [length] = freePathLength(vector, normal, cosI, absorptionCoefficient)
        length = -(1/absorptionCoefficient)*log(rand)*cosI;
    end

    function [reflectedVector] = reflect(vector, normal, cosI)
        % Calculates the reflected vector using angle of incidence
        % equals angle of reflectance
        reflectedVector = vector - normal * 2 * (-cosI);
    end

    function [refracted] = refract( vector, normal, cosI, n1, n2)
        %Computes refracted vector using Snell's law and Heckbert's method
        n = n1/n2;
        %Note: not accounting for total internal reflection on purpose
        %(handled in fresnel)
        refracted = n*vector + (n*cosI - sqrt(1 - n^2*(1-cosI^2)))*normal;
    end

    function [R] = fresnelCoefficient(vector, normal, cosI, n1, n2 )
        %Computed fresnel equation for unpolarized light 
        sinISquared = 1 - cosI^2;
        nSquared = (n1/n2)^2;
        rootTerm = 1-nSquared *sinISquared;
        rootedTerm = sqrt(rootTerm);

        rS = ((n1*cosI - n2*rootedTerm) / (n1*cosI + n2*rootedTerm))^2;
        rP = ((n1*rootedTerm - n2*cosI)/ (n1*rootedTerm +n2*cosI))^2;

        R = (rS + rP) / 2;
        R(rootTerm < 0 ) = 1; %Total internal reflection
    end

    function [ perturbed ] = brakkeScattering(vector, delta)
        % Perturbs a vector according to exponentiated cosine distrubtion,
        % as described by Brakke et al. (1989)
        %[theta, phi, r] = cart2sph(vector(1), vector(2), vector(3)); --
        % WTF matlab
        perturbed = -vector;
        %Rejection sample, perturbing shouldn't flip reflect/refract
        while(sign(perturbed(3)) ~= sign(vector(3)))
            theta = acos(vector(3));
            phi   = atan2(vector(2), vector(1));
            alpha = acos(rand()^(1/(delta + 1)));
            beta = 2*pi*rand();
            perturbed = zeros(1,3);

            newTheta = alpha + theta;
            newPhi   = phi + beta;
            st = sin(newTheta);
            sp = sin(newPhi);
            ct = cos(newTheta);
            cp = cos(newPhi);

            perturbed(1) = st .* cp;
            perturbed(2) = st .* sp;
            perturbed(3) = ct;
        end

        %[x,y,z] = sph2cart(theta + beta, phi+alpha,1);
        %perturbed(1) = x;
        %perturbed(2) = y;
        %perturbed(3) = z;
    end

end

