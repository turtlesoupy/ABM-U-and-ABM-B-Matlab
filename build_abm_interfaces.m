function [interfaces] = build_abm_interfaces(sample, wavelength, bifacial)
    l= load('tissue_lookups.mat');

    interpolationMethod = 'linear';

    dryMatterConcentration = sample.dryBulkDensity / ...
        (1 - sample.airVolumeFraction);
    
    proteinConcentration     = dryMatterConcentration * ...
        sample.proteinFraction;
    celluloseConcentration   = dryMatterConcentration * ...
        sample.celluloseFraction;
    linginConcentration      = dryMatterConcentration * ...
        sample.ligninFraction;
    
    proteinSAC = interp1(l.proteinSACWavelengthLookup, l.proteinSACLookup, wavelength, ...
        interpolationMethod, 'extrap');
    
    celluloseLinginSAC = interp1(l.celluloseLinginSACWavelengthLookup, l.celluloseLinginSACLookup, ...
        wavelength, interpolationMethod, 'extrap');
    
    proteinAbsorptionCoefficient   = proteinConcentration * proteinSAC;
    celluloseAbsorptionCoefficient = celluloseConcentration * celluloseLinginSAC;
    linginAbsorptionCoefficient    = linginConcentration * celluloseLinginSAC;
    
    chlorophyllAbsorption = ...
        (sample.chlorophyllAConcentration + ...
         sample.chlorophyllBConcentration) * ...
         interp1(l.chlorophyllSACWavelengthLookup, l.chlorophyllSACLookup, wavelength, ...
              interpolationMethod, 'extrap');
     
    cartenoidsAbsorption = sample.carotenoidMesophyllConcentration * ...
         interp1(l.caretonoidSACWavelengthLookup, l.caretonoidSACLookup, wavelength, ...
            interpolationMethod, 'extrap');
    
    waterAbsorption = interp1(l.waterSACWavelengthLookup, l.waterSACLookup, wavelength, ...
            interpolationMethod, 'extrap');
      
    mesophyllAbsorption = chlorophyllAbsorption + cartenoidsAbsorption + ...
        proteinAbsorptionCoefficient + celluloseAbsorptionCoefficient + ...
        linginAbsorptionCoefficient + waterAbsorption;
     
    refractiveIndexCuticle =  ...
        interp1(l.refractiveIndexCuticleWavelengthLookup, l.refractiveIndexCuticleLookup, ...
            wavelength,	interpolationMethod, 'extrap');
        
    refractiveIndexWater = ...
         interp1(l.refractiveIndexWaterWavelengthLookup, l.refractiveIndexWaterLookup, ...
            wavelength,	interpolationMethod, 'extrap');
         
    
    refractiveIndexAir = 1;
    refractiveIndexMesophyll = 1.415;
    refractiveIndexAntidermalWall = (1 - sample.antidermalScattererFraction) * ...
        refractiveIndexWater + 1.535 * sample.antidermalScattererFraction;
    
    if bifacial
        interfaces = ABMB_interfaces();
    else
        interfaces = ABMU_interfaces();
    end

    function [interfaces] = ABMU_interfaces()
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
        epidermisMesophyll.perturbanceDownBottom  = sample.spongyCellCapsAspectRatio; 
        epidermisMesophyll.perturbanceUpBottom    = epidermisMesophyll.perturbanceDownBottom;

        mesophyllAir = interfaceStruct;
        mesophyllAir.name = 'Mesophyll<->Air';
        mesophyllAir.n1 = refractiveIndexMesophyll;
        mesophyllAir.n2 = refractiveIndexAir;
        mesophyllAir.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        mesophyllAir.perturbanceUpTop      = mesophyllAir.perturbanceDownTop;
        mesophyllAir.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
        mesophyllAir.perturbanceUpBottom   = mesophyllAir.perturbanceDownBottom;
        mesophyllAir.splitThicknessIndex   = 4;
        mesophyllAir.thickness             = 0.8 * sample.wholeLeafThickness; %Unifacial ratio
        mesophyllAir.absorptionCoefficient  = mesophyllAbsorption;
        
        airMesophyll = interfaceStruct;
        airMesophyll.name = 'Air<->Mesophyll';
        airMesophyll.n1 = refractiveIndexAir;
        airMesophyll.n2 = refractiveIndexMesophyll;
        airMesophyll.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        airMesophyll.perturbanceUpTop      = airMesophyll.perturbanceDownTop;
        airMesophyll.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
        airMesophyll.perturbanceUpBottom   = airMesophyll.perturbanceDownBottom;
        airMesophyll.thickness             = 0.8 * sample.wholeLeafThickness; %Unifacial ratio
        airMesophyll.absorptionCoefficient  = mesophyllAbsorption;

        mesophyllEpidermis = interfaceStruct;
        mesophyllEpidermis.name = 'Mesophyll<->Abaxial Epidermis';
        mesophyllEpidermis.n1   = refractiveIndexMesophyll;
        mesophyllEpidermis.n2   = refractiveIndexCuticle;
        mesophyllEpidermis.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        mesophyllEpidermis.perturbanceUpTop      = mesophyllEpidermis.perturbanceDownTop;
        mesophyllEpidermis.perturbanceDownBottom = sample.epidermisCellCapsAspectRatio;
        mesophyllEpidermis.perturbanceUpBottom   = mesophyllEpidermis.perturbanceDownBottom;

        epidermisAir = interfaceStruct;
        epidermisAir.name = 'Abaxial Epidermis<->Air';
        epidermisAir.n1   = refractiveIndexCuticle;
        epidermisAir.n2   = refractiveIndexAir;
        epidermisAir.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        epidermisAir.perturbanceUpTop      = epidermisAir.perturbanceDownTop;
        epidermisAir.perturbanceUpBottom   = sample.epidermisCellCapsAspectRatio;

        interfaces = repmat(abm_interface(), 1, 6);
        interfaces(1) = airCuticle;
        interfaces(2) = epidermisMesophyll;
        interfaces(3) = mesophyllAir;
        interfaces(4) = airMesophyll;
        interfaces(5) = mesophyllEpidermis;
        interfaces(6) = epidermisAir;
    end

    function [interfaces] = ABMB_interfaces()
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
end

