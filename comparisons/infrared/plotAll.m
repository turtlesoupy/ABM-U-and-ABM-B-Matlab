
% program to plot lopex reflectance and transmittance values
% as well as modeled reflectance and transmittance values (in the visible range) 

% test version to resume testing on 18/03/2010

close all
clear all


oR = load('opex0141');
oT = load('opex0142');
oW = load('opex.wvl');

mR = load('MyR141IRx_DataF.txt');
mT= load('MyT142IRx_DataF.txt');

bR = load('R141IRx.txt');
bT = load('T142IRx.txt');

w = 750:5:2500;


rmsRError = sqrt(mse(bR - mR))
rmsTErorr = sqrt(mse(bT - mT))

wavelengthsIndex = 1:length(w);
matchup = 1 + (wavelengthsIndex - 1) * 5 + 350;

length(mR)
length(matchup)

rmsMRError = sqrt(mse(mR - oR(matchup)))
rmsMTErorr = sqrt(mse(mT - oT(matchup)))

rmsBRError = sqrt(mse(bR - oR(matchup)))
rmsBTErorr = sqrt(mse(bT - oT(matchup)))



subplot(3,2,1)
plot(w,bR*100,'r--',w,mR*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Modelled (Baranoski)','Modelled (Dimson)','Location','NorthEast');
title('(g)','fontsize',12);
axis([750 2500 0 60]);


subplot(3,2,2)
plot(w,bT*100,'r--',w,mT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Modelled (Baranoski)','Modelled (Dimson)','Location','NorthEast');
title('(h)','fontsize',12);
axis([750 2500 0 60]);


subplot(3,2,3)
plot(oW,oR*100,'b-',w,mR*100,'r--','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Measured','Modelled (Dimson)','Location','NorthEast');
title('(g)','fontsize',12);
axis([750 2500 0 60]);


subplot(3,2,4)
plot(oW,oT*100,'b-',w,mT*100,'r--','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Measured','Modelled (Dimson)','Location','NorthEast');
title('(h)','fontsize',12);
axis([750 2500 0 60]);

subplot(3,2,5)
plot(oW,oR*100,'b-',w,bR*100,'r--','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Measured','Modelled (Baranoski)','Location','NorthEast');
title('(g)','fontsize',12);
axis([750 2500 0 60]);


subplot(3,2,6)
plot(oW,oT*100,'b-',w,bT*100,'r--','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Measured','Modelled (Baranoski)','Location','NorthEast');
title('(h)','fontsize',12);
axis([750 2500 0 60]);

pause
%print -depsc fig1-B2.eps
close








