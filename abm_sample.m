function [abm_sample] = abm_sample()
    abm_sample = struct(...
        'wholeLeafThickness',               2.04e-4, ...
        'cuticleUndulationsAspectRatio',    5.0, ...
        'epidermisCellCapsAspectRatio',     5.0, ...
        'spongyCellCapsAspectRatio',        5.0, ...
        'palisadeCellCapsAspectRatio',      1.0, ...
        'proteinConcentration',             0.0, ...
        'celluloseConcentration',           0.0, ...
        'linginConcentration',              0.0, ...
        'chlorophyllAConcentration',        3.978, ...
        'chlorophyllBConcentration',        1.161, ...
        'carotenoidMesophyllConcentration', 1.32, ...
        'mesophyllFraction',                0.8,  ...
        'bifacial',                         0 ...
    );
end

