function [ sample ] = lopex_0141_0142()
    sample = abm_sample();
    sample.wholeLeafThickness        = 2.04e-4;
    
    %in g/cm^3
    sample.linginConcentration       = 0.059245619;
    sample.celluloseConcentration    = 0.0;
    sample.proteinConcentration      = 0.05308714;
    sample.chlorophyllAConcentration = 0.002895146;
    sample.chlorophyllBConcentration = 0.00079866;
    sample.carotenoidConcentration   = 0.000658895;
    
    sample.cuticleUndulationsAspectRatio = 10.0;
    sample.epidermisCellCapsAspectRatio = 5.0;
    sample.spongyCellCapsAspectRatio = 5.0;
    sample.mesophyllFraction         = 0.8;
end

