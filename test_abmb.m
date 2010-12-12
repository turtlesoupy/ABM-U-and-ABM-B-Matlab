sample = abm_sample();

sample.wholeLeafThickness = 2 * 0.83e-4;
%sample.wholeLeafThickness = 1.72;

dryWeight = 1.19e-5;
freshWeight = 4.94e-5;
area = 4.1e-4;

dryBulk = dryWeight / ((1-0.31)*(sample.wholeLeafThickness  * area));
%dryBulk = dryWeight / ((0.5)*(sample.wholeLeafThickness  * area));
freshBulk = freshWeight / ((0.5)*(sample.wholeLeafThickness * area));

sample.linginConcentration       = 0.0424*dryBulk;
sample.celluloseConcentration    = 0.1490*dryBulk;
sample.proteinConcentration      = 0.3106*dryBulk;
sample.chlorophyllAConcentration = (2.74 / 1000) * freshBulk;
sample.chlorophyllBConcentration = (0.8 / 1000) * freshBulk;
sample.carotenoidConcentration   = (0.78 / 1000) * freshBulk;
sample.mesophyllFraction         = 0.5;

sample.cuticleUndulationsAspectRatio = 5;
sample.epidermisCellCapsAspectRatio = 5;
sample.spongyCellCapsAspectRatio = 5;
sample.palisadeCellCapsAspectRatio = 1;
sample.bifacial = 1;


nSamples = 100000;
incidentAngle = 172 * pi / 180;
wavelengths = 400e-9:5e-9:2500e-9;

nWavelengths = length(wavelengths);
reflectance = zeros(1, nWavelengths);
transmittance = zeros(1, nWavelengths);
absorptance = zeros(1, nWavelengths);

for n = 1:length(wavelengths)
    tic
    fprintf('Wavelength %d\n', wavelengths(n))
    
    interfaces = build_interfaces(sample, wavelengths(n));
    %sac(n) = interfaces(2).absorptionBelow;
    %for i = 1:length(interfaces)
    %    interfaces(i)
    %end
    [r, t, a] = ABM(0, incidentAngle, interfaces,nSamples);
    reflectance(n) = r;
    transmittance(n) = t;
    absorptance(n) = a;
    toc
end 

save('abm_test')