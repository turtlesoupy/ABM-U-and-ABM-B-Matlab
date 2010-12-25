addpath samples src;
sample = lopex_0141_0142();

samplesPerWavelength  = 1000;
wavelengths           = 400e-9:5e-9:2500e-9;
polarAngle            = 8 * pi / 180;
azimuthalAngle        = 0;

[reflectance, transmittance, absorptance] = ...
    abmu(sample, samplesPerWavelength, wavelengths, azimuthalAngle, polarAngle,0);

subplot(2,1,1)
plot(wavelengths * 10e8, reflectance*100,'r--','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
title('ABM-U','fontsize',12);
axis([400 2500 0 60]);

subplot(2,1,2)
plot(wavelengths * 10e8, transmittance*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
title('ABM-U','fontsize',12);
axis([400 2500 0 60]);