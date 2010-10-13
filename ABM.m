function [reflectance, transmittance, absorptance] = ABM(azimuthalI, polarI, interfaceArray, nSamples)
    stepFunction  = @step;     
    [x,y,z] = sph2cart(azimuthalI, -polarI + pi/2, 1); %Match canonical
    startDirection = [x,y,z];
    
    if startDirection(3) < 0 
         startState    = 1;
         reflectedState = 0;
         transmittedState = length(interfaceArray) + 1;
    else
        startState    = length(interfaceArray);  
        reflectedState = length(interfaceArray) + 1;
        transmittedState = 0;
    end
    
    absorbedState   = length(interfaceArray) + 2;
    endStates       = zeros(1, nSamples);
    startSplits     = find([interfaceArray.splitThicknessIndex]);
    endSplits       = [interfaceArray(startSplits).splitThicknessIndex];
    numSplits       = length(startSplits);
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
        while(state ~= reflectedState && state ~= transmittedState && state ~= absorbedState)
            [state, direction] = stepFunction(state, direction, thicknesses);
        end
        
        endStates(i) = state;
    end
    
    reflectance   = sum(endStates == reflectedState) / nSamples;
    transmittance = sum(endStates == transmittedState) / nSamples;
    absorptance   = sum(endStates == absorbedState) / nSamples;
    
    function [outState, outVector] = step(state, vector, interfaceThicknesses)
        interface = interfaceArray(state);
        thickness = interfaceThicknesses(state);
      %  state
        if vector(3) < 0 
            normal = [0,0,1];
            n1 = interface.n1;
            n2 = interface.n2;
            perturbanceReflect = interface.perturbanceDownTop;
            perturbanceRefract = interface.perturbanceDownBottom;
            refractState = state + 1;
            reflectState = state - 1;
        else
            normal = [0,0,-1];
            n1 = interface.n2;
            n2 = interface.n1;
            perturbanceReflect = interface.perturbanceUpBottom;
            perturbanceRefract = interface.perturbanceUpTop;
            reflectState = state + 1;
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
                outState = reflectState;
                outVector = reflect(vector, normal, normalAngle);
                if ~isinf(perturbanceReflect)
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

    function [cart] = toCart(azimuthal, polar, basis)
        sa = sin(azimuthal);
        sp = sin(polar);
        ca = cos(azimuthal);
        cp = cos(polar);
        cart = zeros(1,3);
        cart(1) = sp .* ca;
        cart(2) = sp .* sa;
        cart(3) = cp;
    end

    function [u,v,w] = basis(vector)
        u = vector / norm(vector);
        u = u/norm(u);
        w = cross(u, perpendicular(vector));
        w = w/norm(w);
        v = cross(w,u);     
        v = v / norm(v);
    end

    function [perp] = perpendicular(vector)
        perp = [vector(2), -vector(1), 0];
        if abs(dot(vector, perp)) < 0.01
            perp = [0, vector(3), -vector(2)];
        end
    end

    function [ perturbed ] = brakkeScattering(vector, delta)
        % Perturbs a vector according to exponentiated cosine distrubtion,
        % as described by Brakke et al. (1989)

        if abs(vector(1)) < abs(vector(2)) && abs(vector(1)) < abs(vector(3))
            perp = [0, -vector(3), vector(2)];
        elseif abs(vector(2)) < abs(vector(3))
            perp = [vector(3), 0, -vector(1)];
        else
            perp = [-vector(2), vector(1), 0];
        end

        w = vector / norm(vector);
        u = perp / norm(perp);
        % u = u - w * dot(w,u)
        % u = u / norm(u);
        v = cross(w,u);
        perturbed = -vector;
        i = 0;
        while(sign(perturbed(3)) ~= sign(vector(3)))
           polar = acos(rand^(1/(delta + 1)));
           azimuthal  = 2*pi*rand;
           sp = sin(polar);
           sa = sin(azimuthal);
           cp = cos(polar);
           ca = cos(azimuthal);

           perturbed = u * (sp * ca) + (v * sp * sa) + (w * cp);  
           i = i + 1;
           if i >= 100
               fprintf('Broke infinite loop for perturbing\n');
               perturbed = vector;
               return;
           end
        end
    end

end

