function [ sample ] = lopex_0219_0220()
    sample = abm_sample();
    sample.wholeLeafThickness            = 1.66e-4;
    
    sample.linginConcentration       = 10.7441;
    sample.celluloseConcentration    = 37.7565;
    sample.proteinConcentration      = 78.7059;
    sample.chlorophyllAConcentration = 3.9775;
    sample.chlorophyllBConcentration = 1.1613;
    sample.carotenoidConcentration   = 1.1323;

    sample.cuticleUndulationsAspectRatio = 5.0;
    sample.epidermisCellCapsAspectRatio  = 5.0;
    sample.spongyCellCapsAspectRatio     = 5.0;
    sample.palisadeCellCapsAspectRatio   = 1.0;
    sample.mesophyllFraction         = 0.5;
end

