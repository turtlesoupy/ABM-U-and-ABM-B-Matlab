%Test ABM
sample = abm_sample();
fw = 6.88e-5;
dw = 1.46e-5;
th = 2.04e-4;
area = 4.1e-4; %* 0.98;
obs = 5.0; %* 1.2;

%Following my TestABM script
%th = th * 0.8;
%fw = fw * 0.75 + 0.25 * dw;

mesophyllFreshBulkDensity = fw  / (4.1e-4 * th * 0.8);

sample.wholeLeafThickness = th;
sample.dryBulkDensity = dw  / (4.1e-4 * th * 0.8);
% non-af sample.chlorophyllAConcentration        = 0.003   * mesophyllFreshBulkDensity;
% non-af sample.chlorophyllBConcentration        = 0.00083 * mesophyllFreshBulkDensity;
% AF sample.chlorophyllAConcentration        = 0.0029   * mesophyllFreshBulkDensity;
% AF sample.chlorophyllBConcentration        = 0.0008 * mesophyllFreshBulkDensity;
sample.chlorophyllAConcentration        = 0.0029   * mesophyllFreshBulkDensity;
sample.chlorophyllBConcentration        = 0.0008 * mesophyllFreshBulkDensity;
sample.carotenoidMesophyllConcentration = 0.00066 * mesophyllFreshBulkDensity * 4;
%non-af sample.proteinFraction = 0.2531;
%non-af sample.ligninFraction = 0.2399;
%AF sample.proteinFraction = 0.2655;
%AF sample.ligninFraction = 0.2963;
sample.proteinFraction = 0.2655;
sample.ligninFraction = 0.2963;

sample.cuticleUndulationsAspectRatio = 10;
sample.epidermisCellCapsAspectRatio = 5;
sample.spongyCellCapsAspectRatio = 5;
sample.airVolumeFraction = 0;
sample.spongyCellCapsAspectRatio = obs;
nSamples = 100000;
wavelengths = 695e-9:5e-9:700e-9;
nWavelengths = length(wavelengths);
reflectance = zeros(1, nWavelengths);
transmittance = zeros(1, nWavelengths);
absorptance = zeros(1, nWavelengths);
tic

incidentAngle = 8 * pi / 180;

for n = 1:length(wavelengths)
    tic
    fprintf('Wavelength %d\n', wavelengths(n))
    
    interfaces = build_abm_interfaces(sample, wavelengths(n), 0);
    [r, t, a] = ABM(0, incidentAngle, interfaces,nSamples);
    reflectance(n) = r;
    transmittance(n) = t;
    absorptance(n) = a;
    toc
end
toc

wavelengths
reflectance
transmittance
absorptance