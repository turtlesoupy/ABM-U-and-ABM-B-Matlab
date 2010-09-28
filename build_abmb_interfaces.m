function [interfaces] = build_abmb_interfaces(sample, wavelength)
    load tissue_lookups.mat
    %ABM-B
    interpolationMethod = 'linear';

    dryMatterConcentration = sample.dryBulkDensity / ...
        (1 - sample.airVolumeFraction);
    
    proteinConcentration     = dryMatterConcentration * ...
        sample.proteinFraction;
    celluloseConcentration   = dryMatterConcentration * ...
        sample.celluloseFraction;
    linginConcentration      = dryMatterConcentration * ...
        sample.ligninFraction;
    
    proteinSAC = interp1(proteinSACWavelengthLookup, proteinSACLookup, wavelength, ...
        interpolationMethod, 'extrap');
    
    celluloseLinginSAC = interp1(celluloseLinginSACWavelengthLookup, celluloseLinginSACLookup, ...
        wavelength, interpolationMethod, 'extrap');
    
    proteinAbsorptionCoefficient   = proteinConcentration * proteinSAC;
    celluloseAbsorptionCoefficient = celluloseConcentration * celluloseLinginSAC;
    linginAbsorptionCoefficient    = linginConcentration * celluloseLinginSAC;
    chlorophyllAbsorption = ...
        (sample.chlorophyllAConcentration + ...
         sample.chlorophyllBConcentration) * ...
         interp1(chlorophyllSACWavelengthLookup, chlorophyllSACLookup, wavelength, ...
              interpolationMethod, 'extrap');
     
    cartenoidsAbsorption = sample.carotenoidMesophyllConcentration * ...
         interp1(caretonoidSACWavelengthLookup, caretonoidSACLookup, wavelength, ...
            interpolationMethod, 'extrap');
    
    waterAbsorption = interp1(waterSACWavelengthLookup, waterSACLookup, wavelength, ...
            interpolationMethod, 'extrap');
      
    mesophyllAbsorption = chlorophyllAbsorption + cartenoidsAbsorption + ...
        proteinAbsorptionCoefficient + celluloseAbsorptionCoefficient + ...
        linginAbsorptionCoefficient + waterAbsorption;
     
    refractiveIndexCuticle =  ...
        interp1(refractiveIndexCuticleWavelengthLookup, refractiveIndexCuticleLookup, ...
            wavelength,	interpolationMethod, 'extrap');
        
    refractiveIndexWater = ...
         interp1(refractiveIndexWaterWavelengthLookup, refractiveIndexWaterLookup, ...
            wavelength,	interpolationMethod, 'extrap');
         
    
    refractiveIndexAir = 1;
    refractiveIndexMesophyll = 1.415;
    refractiveIndexAntidermalWall = (1 - sample.antidermalScattererFraction) * ...
        refractiveIndexWater + 1.535 * sample.antidermalScattererFraction;
    
    
    %Interfaces Declaration
    interfaceStruct = abm_interface();

    airCuticle = interfaceStruct;
    airCuticle.name = 'Air<->Adaxial Epidermis';
    airCuticle.n1 = refractiveIndexAir;
    airCuticle.n2 = refractiveIndexCuticle;
    airCuticle.perturbanceDownTop    = sample.cuticleUndulationsAspectRatio;
    airCuticle.perturbanceDownBottom = sample.epidermisCellCapsAspectRatio;
    airCuticle.perturbanceUpBottom   = sample.epidermisCellCapsAspectRatio;
   
    
    epidermisMesophyll = interfaceStruct;
    epidermisMesophyll.name = 'Adaxial Epidermis<->Mesophyll';
    epidermisMesophyll.n1 = refractiveIndexCuticle;
    epidermisMesophyll.n2 = refractiveIndexMesophyll;
    epidermisMesophyll.perturbanceDownTop     = sample.epidermisCellCapsAspectRatio;
    epidermisMesophyll.perturbanceUpTop       = epidermisMesophyll.perturbanceDownTop;
    epidermisMesophyll.perturbanceDownBottom  = sample.palisadeCellCapsAspectRatio;
    epidermisMesophyll.perturbanceUpBottom    = epidermisMesophyll.perturbanceDownBottom;
    epidermisMesophyll.thickness              = 0.5 * sample.wholeLeafThickness; %Bifacial ratio
    epidermisMesophyll.absorptionCoefficient  = mesophyllAbsorption;
    
    mesophyllAir = interfaceStruct;
    mesophyllAir.name = 'Mesophyll<->Air';
    mesophyllAir.n1 = refractiveIndexMesophyll;
    mesophyllAir.n2 = refractiveIndexAir;
    mesophyllAir.perturbanceDownTop    = sample.palisadeCellCapsAspectRatio;
    mesophyllAir.perturbanceUpTop      = mesophyllAir.perturbanceDownTop;
    mesophyllAir.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
    mesophyllAir.perturbanceUpBottom   = mesophyllAir.perturbanceDownBottom;
    
    
    airAntidermalWall = interfaceStruct;
    airAntidermalWall.name = 'Air<->Antidermal Wall';
    airAntidermalWall.n1 = refractiveIndexAir;
    airAntidermalWall.n2 = refractiveIndexAntidermalWall;
    
    antidermalWallCuticle = interfaceStruct;
    antidermalWallCuticle.name = 'Antidermal Wall<->Abaxial Epidermis';
    antidermalWallCuticle.n1 = refractiveIndexAntidermalWall;
    antidermalWallCuticle.n2 = refractiveIndexCuticle;
    antidermalWallCuticle.perturbanceDownBottom = sample.epidermisCellCapsAspectRatio;
    antidermalWallCuticle.perturbanceUpBottom = antidermalWallCuticle.perturbanceDownBottom;
    
    cuticleAir = interfaceStruct;
    cuticleAir.name = 'Abaxial Epidermis<->Air';
    cuticleAir.n1 = refractiveIndexCuticle;
    cuticleAir.n2 = refractiveIndexAir;
    cuticleAir.perturbanceDownTop  = sample.epidermisCellCapsAspectRatio;
    cuticleAir.perturbanceUpTop    = cuticleAir.perturbanceDownTop;
    cuticleAir.perturbanceUpBottom = sample.cuticleUndulationsAspectRatio;
    
    
    interfaces = repmat(airCuticle, 1, 6);
    interfaces(1) = airCuticle;
    interfaces(2) = epidermisMesophyll;
    interfaces(3) = mesophyllAir;
    interfaces(4) = airAntidermalWall;
    interfaces(5) = antidermalWallCuticle;
    interfaces(6) = cuticleAir;
end

