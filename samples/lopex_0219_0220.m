function [ sample ] = lopex_0219_0220()
    sample = abm_sample();
    sample.wholeLeafThickness            = 1.66e-4;
    
    %in g/cm^3
    sample.linginConcentration       = 0.0107441;
    sample.celluloseConcentration    = 0.0377565;
    sample.proteinConcentration      = 0.0787059;
    sample.chlorophyllAConcentration = 0.0039775;
    sample.chlorophyllBConcentration = 0.0011613;
    sample.carotenoidConcentration   = 0.0011323;

    sample.cuticleUndulationsAspectRatio = 5.0;
    sample.epidermisCellCapsAspectRatio  = 5.0;
    sample.spongyCellCapsAspectRatio     = 5.0;
    sample.palisadeCellCapsAspectRatio   = 1.0;
    sample.mesophyllFraction         = 0.5;
end

