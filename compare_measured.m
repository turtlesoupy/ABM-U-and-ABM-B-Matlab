                                                                     
                                                                     
                                                                     
                                             
% program to plot lopex reflectance and transmittance values
% as well as modeled reflectance and transmittance values (in the visible range) 

% test version to resume testing on 18/03/2010

close all
clear all

op = load('opex.wvl');
%op=opex;
op
r = load('opex0141');
%r=opex0141;

t = load('opex0142');
%t=opex0142;

m = load('MyR141AF.txt');
%m=R141AF;

n = load('MyT142AF.txt');
%n=T142af;

w(1)=400;

for i=2:61
    w(i)=w(i-1)+5;
end;

maxe=0;
sume=0;
maxae=0;
sumae=0;
sumae2=0;
for i=1:61
    e=abs(r(1+(i-1)*5)-m(i))/r(1+(i-1)*5);
    ae=abs(r(1+(i-1)*5)-m(i));
    sume=sume+e;
    sumae=sumae+ae;
    sumae2=sumae2+(ae*ae);
    if (e>maxe)
        maxe=e;
        j=i;
        ind=400+(i-1)*5;
    end
    if (ae>maxae)
        maxae=ae;
        ja=i;
        inda=400+(i-1)*5;
    end
end
rmse=sqrt(sumae2/61)
avge=sume/61;
avgae=sumae/61;

maxe=0;
sume=0;
maxae=0;
sumae=0;
sumae2=0;
for i=1:61
    e=abs(t(1+(i-1)*5)-n(i))/t(1+(i-1)*5);
    ae=abs(t(1+(i-1)*5)-n(i));
    sume=sume+e;
    sumae=sumae+ae;
    sumae2=sumae2+(ae*ae);
    if (e>maxe)
        maxe=e;
        j=i;
        ind=400+(i-1)*5;
    end
    if (ae>maxae)
        maxae=ae;
        ja=i;
        inda=400+(i-1)*5;
    end
end
rmse=sqrt(sumae2/61)
avge=sume/61;
avgae=sumae/61;


subplot(3,2,1)
plot(op,r*100,'r--',w,m*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('reflectance (%)','fontsize',12);
legend('measured','modeled','Location','NorthWest');
title('(g)','fontsize',12);
axis([400 700 5 25]);


subplot(3,2,2)
plot(op,t*100,'r--',w,n*100,'b-','linewidth',2);
xlabel('wavelength (nm)','fontsize',12);
ylabel('transmittance (%)','fontsize',12);
legend('measured','modeled','Location','NorthWest');
title('(h)','fontsize',12);
axis([400 700 0 20]);


pause
%print -depsc fig1-B2.eps
close