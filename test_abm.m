sample = abm_sample();

sample.linginConcentration       = 59.245619;
sample.celluloseConcentration    = 0;
sample.proteinConcentration      = 53.08714;
sample.chlorophyllAConcentration = 2.895146;
sample.chlorophyllBConcentration = 0.79866;
sample.carotenoidConcentration   = 0.658895;
sample.mesophyllFraction         = 0.8;

sample.cuticleUndulationsAspectRatio = 10;
sample.epidermisCellCapsAspectRatio = 5;
sample.spongyCellCapsAspectRatio = 5;

nSamples = 10000;
incidentAngle = 8 * pi / 180;
wavelengths = 400e-9:5e-9:2500e-9;

nWavelengths = length(wavelengths);
reflectance = zeros(1, nWavelengths);
transmittance = zeros(1, nWavelengths);
absorptance = zeros(1, nWavelengths);

for n = 1:length(wavelengths)
    tic
    fprintf('Wavelength %d\n', wavelengths(n))
    
    interfaces = build_interfaces(sample, wavelengths(n));
   % for i = 1:length(interfaces)
   %     interfaces(i)
   % end
    [r, t, a] = ABM(0, incidentAngle, interfaces,nSamples);
    reflectance(n) = r;
    transmittance(n) = t;
    absorptance(n) = a;
    toc
end 
