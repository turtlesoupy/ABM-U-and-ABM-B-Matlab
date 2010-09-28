%Test ABM
sample = abm_sample();
fw = 6.88e-5;
dw = 1.46e-5;
th = 2.04e-4;

area = 4.1e-4 * 0.98;
obs = 5.0 * 1.2;

%Following my TestABM script
th = th * 0.8;
fw = fw * 0.75 + 0.25 * dw;

mesophyllFreshBulkDensity = fw  * 4.1e-4 * th * 0.8;

sample.wholeLeafThickness = th;
sample.dryBulkDensity = dw  / (4.1e-4 * th * 0.8);
sample.chlorophyllAConcentration        = 0.003   * mesophyllFreshBulkDensity;
sample.chlorophyllBConcentration        = 0.00083 * mesophyllFreshBulkDensity;
sample.carotenoidMesophyllConcentration = 0.00066 * mesophyllFreshBulkDensity * 4;
sample.proteinFraction = 0.2531;
sample.ligninFraction = 0.2399;
sample.cuticleUndulationsAspectRatio = 10;
sample.epidermisCellCapsAspectRatio = 5;
sample.spongyCellCapsAspectRatio = 5;
sample.airVolumeFraction = 0;
sample.spongyCellCapsAspectRatio = obs;

interfaces = build_abmb_interfaces(sample, 400e-9)
nSamples = 10000
%'New code'
%tic
%[r, t, a] = ABM(0, -30 * pi/180, interfaces,nSamples)
%toc
'Old Code'
tic
[r, t, a] = ABM_pre_vectorize(0, -30 * pi/180, interfaces,nSamples)
toc