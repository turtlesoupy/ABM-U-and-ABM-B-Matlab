
% program to plot lopex reflectance and transmittance values
% as well as modeled reflectance and transmittance values (in the visible range) 

% test version to resume testing on 18/03/2010

close all
clear all

op = load('opex.wvl');
r = load('opex0141');
t= load('opex0142');

mR = load('MyR141AF-2.txt');
mT= load('MyT142AF-2.txt');

bR = load('R141AF.txt');
bT = load('T142AF.txt');

sqrt(mse(bR - mR))

w = 400:5:700;

wavelengthsIndex = 1:length(w);
matchup = 1 + (wavelengthsIndex - 1) * 5;

rmsBReflectance   = sqrt(mse(r(matchup) - bR))
rmsBTransmittance = sqrt(mse(t(matchup) - bT))
rmsMReflectance   = sqrt(mse(r(matchup) - mR))
rmsMTransmittance = sqrt(mse(t(matchup) - mT))

subplot(3,2,1)
plot(op,r*100,'r--',w,bR*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Measured','Modeled (Baranoski)','Location','NorthWest');
title('(g)','fontsize',12);
axis([400 700 5 25]);


subplot(3,2,2)
plot(op,t*100,'r--',w,bT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Measured','Modeled (Baranoski)','Location','NorthWest');
title('(h)','fontsize',12);
axis([400 700 0 20]);


subplot(3,2,3)
plot(op,r*100,'r--',w,mR*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('Measured','Modeled (Dimson)','Location','NorthWest');
title('(g)','fontsize',12);
axis([400 700 5 25]);

subplot(3,2,4)
plot(op,t*100,'r--',w,mT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Measured','Modeled (Dimson)','Location','NorthWest');
title('(h)','fontsize',12);
axis([400 700 0 20]);

subplot(3,2,5)
plot(w,bT*100,'r--',w,mT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Measured (Baranoski)','Modeled (Dimson)','Location','NorthWest');
title('(h)','fontsize',12);
axis([400 700 0 20]);

subplot(3,2,6)
plot(w,bT*100,'r--',w,mT*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('Modelled (Baranoski)','Modeled (Dimson)','Location','NorthWest');
title('(h)','fontsize',12);
axis([400 700 0 20]);

pause
%print -depsc fig1-B2.eps
close








