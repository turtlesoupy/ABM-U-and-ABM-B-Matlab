
% program to plot lopex reflectance and transmittance values
% as well as modeled reflectance and transmittance values (in the visible range) 

% test version to resume testing on 18/03/2010

close all
clear all

mR = load('MyR141IRx_DataF.txt');
mT= load('MyT142IRx_DataF.txt');

bR = load('R141IRx.txt');
bT = load('T142IRx.txt');

w = 750:5:2500;

rmsRError = sqrt(mse(bR - mR))
rmsTErorr = sqrt(mse(bT - mT))

wavelengthsIndex = 1:length(w);
matchup = 1 + (wavelengthsIndex - 1) * 5 + 350;

figure
plot(w,bR*100,'r--',w,mR*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Modelled (Baranoski)','Modelled (Dimson)','Location','NorthEast');
title('Reflectance (%) comparisons (100,000 rays, incidence [0.1392, 0, 0.9903], normal [0,0,1], data set F))','fontsize',12);
axis([750 2500 0 60]);

figure
plot(w,bT*100,'r--',w,mT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Modelled (Baranoski)','Modelled (Dimson)','Location','NorthEast');
title('Transmittance (%) comparisons (100,000 rays, incidence [0.1392, 0, 0.9903], normal top [0,0,1], data set F))','fontsize',12);
axis([750 2500 0 60]);