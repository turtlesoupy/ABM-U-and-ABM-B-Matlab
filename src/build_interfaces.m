function [interfaces] = build_interfaces(sample, wavelength, bifacial)
    [PATHSTR, NAME, EXT] = fileparts(which('build_interfaces.m'));
    DATA_DIR           = [PATHSTR,  '/../data/'];
    CAROTENOIDS_FILE   = [DATA_DIR, 'caro-PAS-400-2500.txt'];
    CELLULOSE_FILE     = [DATA_DIR, 'cellulose400-2500.txt'];
    CHLOROPHYLL_FILE   = [DATA_DIR, 'chloAB-DFA-400-2500.txt'];
    PROTEIN_FILE       = [DATA_DIR, 'protein400-2500.txt'];
    MESOPHYLL_RI_FILE  = [DATA_DIR, 'rmH400-2500.txt'];
    WATER_SAC_FILE     = [DATA_DIR, 'sacwH400-2500.txt'];
    CUTICLE_RI_FILE    = [DATA_DIR, 'rcH400-2500.txt'];
    ANTIDERMAL_RI_FILE = [DATA_DIR, 'raH400-2500.txt'];

    DATA_WAVELENGTHS = 400e-9:5e-9:2500e-9;
    
    CAROTENOIDS_ABSORPTION_LOOKUP = load(CAROTENOIDS_FILE)' ./ 10;
    CELLULOSE_ABSORPTION_LOOKUP = load(CELLULOSE_FILE)' ./ 10;
    CHLOROPHYLL_ABSORPTION_LOOKUP = load(CHLOROPHYLL_FILE)' ./ 10;
    PROTEIN_ABSORPTION_LOOKUP = load(PROTEIN_FILE)' ./ 10;
    WATER_SAC_LOOKUP = load(WATER_SAC_FILE)' .* 100;
    MESOPHYLL_RI_LOOKUP = load(MESOPHYLL_RI_FILE);
    CUTICLE_RI_LOOKUP = load(CUTICLE_RI_FILE)';
    ANTIDERMAL_RI_LOOKUP = load(ANTIDERMAL_RI_FILE)';

    interpFunc = @(lookup) interp1(DATA_WAVELENGTHS, lookup, wavelength, ... 
                    'linear', 'extrap');

    proteinAbsorptionCoefficient     = sample.proteinConcentration * interpFunc(PROTEIN_ABSORPTION_LOOKUP);
    chlorophyllAbsorptionCoefficient = (sample.chlorophyllAConcentration + sample.chlorophyllBConcentration) ...
                                        * interpFunc(CHLOROPHYLL_ABSORPTION_LOOKUP);
    carotenoidAbsorptionCoefficient  = sample.carotenoidConcentration * interpFunc(CAROTENOIDS_ABSORPTION_LOOKUP);
    celluloseAbsorptionCoefficient   = sample.celluloseConcentration  * interpFunc(CELLULOSE_ABSORPTION_LOOKUP);
    linginAbsorptionCoefficient      = sample.linginConcentration * interpFunc(CELLULOSE_ABSORPTION_LOOKUP);
    
    waterAbsorptionCoefficient = interpFunc(WATER_SAC_LOOKUP);

    mesophyllAbsorption = chlorophyllAbsorptionCoefficient + ...
                          carotenoidAbsorptionCoefficient + ...
                          celluloseAbsorptionCoefficient + ...
                          proteinAbsorptionCoefficient  + ...
                          linginAbsorptionCoefficient + ...
                          waterAbsorptionCoefficient;
    
    refractiveIndexCuticle        = interpFunc(CUTICLE_RI_LOOKUP);
    refractiveIndexMesophyll      = interpFunc(MESOPHYLL_RI_LOOKUP);
    refractiveIndexAntidermalWall = interpFunc(ANTIDERMAL_RI_LOOKUP);
    refractiveIndexAir = 1;

    mesophyllThickness = sample.mesophyllFraction * sample.wholeLeafThickness;
    
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
        epidermisMesophyll.thicknessBelow         = mesophyllThickness;
        epidermisMesophyll.absorptionBelow        = mesophyllAbsorption;        
        epidermisMesophyll.splitBelow             = 1;

        mesophyllAir = interfaceStruct;
        mesophyllAir.name = 'Mesophyll<->Air';
        mesophyllAir.n1 = refractiveIndexMesophyll;
        mesophyllAir.n2 = refractiveIndexAir;
        mesophyllAir.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        mesophyllAir.perturbanceUpTop      = mesophyllAir.perturbanceDownTop;
        mesophyllAir.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
        mesophyllAir.perturbanceUpBottom   = mesophyllAir.perturbanceDownBottom;
        mesophyllAir.thicknessAbove        = mesophyllThickness;
        mesophyllAir.absorptionAbove       = mesophyllAbsorption;
        mesophyllAir.splitAbove            = 1;
        
        airMesophyll = interfaceStruct;
        airMesophyll.name = 'Air<->Mesophyll';
        airMesophyll.n1 = refractiveIndexAir;
        airMesophyll.n2 = refractiveIndexMesophyll;
        airMesophyll.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        airMesophyll.perturbanceUpTop      = airMesophyll.perturbanceDownTop;
        airMesophyll.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
        airMesophyll.perturbanceUpBottom   = airMesophyll.perturbanceDownBottom;
        airMesophyll.thicknessBelow        = mesophyllThickness;
        airMesophyll.absorptionBelow       = mesophyllAbsorption;
        airMesophyll.splitBelow            = 2;

        mesophyllEpidermis = interfaceStruct;
        mesophyllEpidermis.name = 'Mesophyll<->Abaxial Epidermis';
        mesophyllEpidermis.n1   = refractiveIndexMesophyll;
        mesophyllEpidermis.n2   = refractiveIndexCuticle;
        mesophyllEpidermis.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        mesophyllEpidermis.perturbanceUpTop      = mesophyllEpidermis.perturbanceDownTop;
        mesophyllEpidermis.perturbanceDownBottom = sample.epidermisCellCapsAspectRatio;
        mesophyllEpidermis.perturbanceUpBottom   = mesophyllEpidermis.perturbanceDownBottom;
        mesophyllEpidermis.thicknessAbove        = mesophyllThickness;
        mesophyllEpidermis.absorptionAbove       = mesophyllAbsorption;
        mesophyllEpidermis.splitAbove            = 2;

        epidermisAir = interfaceStruct;
        epidermisAir.name = 'Abaxial Epidermis<->Air';
        epidermisAir.n1   = refractiveIndexCuticle;
        epidermisAir.n2   = refractiveIndexAir;
        epidermisAir.perturbanceDownTop    = sample.spongyCellCapsAspectRatio;
        epidermisAir.perturbanceUpTop      = epidermisAir.perturbanceDownTop;
        epidermisAir.perturbanceUpBottom   = sample.cuticleUndulationsAspectRatio;
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
        epidermisMesophyll.thicknessBelow         = mesophyllThickness;
        epidermisMesophyll.absorptionBelow        = mesophyllAbsorption;
        
        mesophyllAir = interfaceStruct;
        mesophyllAir.name = 'Mesophyll<->Air';
        mesophyllAir.n1 = refractiveIndexMesophyll;
        mesophyllAir.n2 = refractiveIndexAir;
        mesophyllAir.perturbanceDownTop    = sample.palisadeCellCapsAspectRatio;
        mesophyllAir.perturbanceUpTop      = mesophyllAir.perturbanceDownTop;
        mesophyllAir.perturbanceDownBottom = sample.spongyCellCapsAspectRatio;
        mesophyllAir.perturbanceUpBottom   = mesophyllAir.perturbanceDownBottom;
        mesophyllAir.absorptionAbove       = mesophyllAbsorption;
        mesophyllAir.thicknessAbove        = mesophyllThickness;
        
        airAntidermalWall = interfaceStruct;
        airAntidermalWall.name = 'Air<->Antidermal Wall';
        airAntidermalWall.n1 = refractiveIndexAir;
        airAntidermalWall.n2 = refractiveIndexAntidermalWall;
        
        antidermalWallCuticle = interfaceStruct;
        antidermalWallCuticle.name = 'Antidermal Wall<->Abaxial Epidermis';
        antidermalWallCuticle.n1 = refractiveIndexAntidermalWall;
        antidermalWallCuticle.n2 = refractiveIndexCuticle;
        antidermalWallCuticle.perturbanceDownBottom = sample.epidermisCellCapsAspectRatio;
        antidermalWallCuticle.perturbanceUpBottom   = antidermalWallCuticle.perturbanceDownBottom;
        
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

