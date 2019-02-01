%%
%  Comparison between non-uniformly spaced array and uniformly spaced array 
%
%  Version 1
%  Zhengyu Peng
%
%%
clear;

wavelength=1;
k=2*pi/wavelength;   % wave number
angleStep=0.05;
theta=0:angleStep:360;
elementNumber=8;
thetaM=50.5;      % location of the main lobe

%% non-uniform array
d1=[0,1.17611290322581,1.68447580645161,2.35025000000000,3.40133064516129,4.02362903225806,5.80238709677419,6.86625806451613];
%w1=[0.104858767920655 + 0.00638377758385419i,-0.155896761387194 - 0.112022190087344i,-0.125230197364286 + 0.103746484301190i,0.0204564334818207 + 0.143544730805173i,0.0510063282098636 - 0.0811331430465989i,-0.115769238179036 - 0.125985531777502i,0.0591959384958559 + 0.0121478103967238i,-0.0585918747674234 - 0.0506620304236720i].';
w1=[0.0641592867347555 - 0.00949813218645218i,0.0213414496230591 + 0.171004670288957i,0.126305559844362 - 0.125331188457008i,-0.140077853533412 - 0.00898555473933959i,0.0290687944868934 - 0.0721485142088897i,-0.148705460650750 + 0.0589554743071804i,-0.0422836884253894 + 0.102505871797712i,-0.0955462768022038 - 0.0555310933213185i].';
%% normal dense array
d2=[0,0.5*ones(1,elementNumber-1)];
for nn=2:length(d2)
    d2(nn)=d2(nn-1)+d2(nn);
end
%w2=[0.206663964404846 - 0.0112377446515177i,0.0391623588928282 - 0.0795400851563046i,-0.0563062985262089 - 0.0820458561893604i,-0.104843473925938 + 0.00942258776671084i,-0.0431344849835628 + 0.0960226816449238i,0.0591681008950777 + 0.0800066097771803i,0.0879672179410876 - 0.0110491659607346i,0.0780914556024583 - 0.191671608879013i].';

w2=chebwin(elementNumber,9.5).*(exp(-1i*k*d2*cosd(thetaM)).');

%% sparse array
d3=[0,d1(elementNumber)/(elementNumber-1)*ones(1,elementNumber-1)];
for nn=2:length(d3)
    d3(nn)=d3(nn-1)+d3(nn);
end
w3=exp(-1i*k*d3*cosd(thetaM)).';

%% 26-element normal dense array
eleNum=12;
d4=[0,0.5*ones(1,eleNum-1)];
for nn=2:length(d4)
    d4(nn)=d4(nn-1)+d4(nn);
end
%w4=[0.151762138140478 - 0.0102832612349678i;-0.0148958031630125 - 0.0172358745075817i;-0.00512570908792318 + 0.0237363670087733i;0.0233375004828243 - 0.0107751751748616i;-0.0238137186663838 - 0.0127944992034719i;0.00428503645225367 + 0.0279267418659784i;0.0205014091227646 - 0.0210035366464110i;-0.0300694264463313 - 0.00387332943146368i;0.0153922908579175 + 0.0270630607817905i;0.0126321128110758 - 0.0291921038237749i;-0.0314120648601946 + 0.00757882283295584i;0.0252090020191759 + 0.0207567876802740i;0.00142889348747491 - 0.0327916992328431i;-0.0270580598489617 + 0.0185795229321596i;0.0309579586761419 + 0.0103895530323716i;-0.0102082265061976 - 0.0306585793198716i;-0.0180712108839995 + 0.0261759534099189i;0.0311078847753170 - 0.00127725546942699i;-0.0192111629602375 - 0.0234543024925927i;-0.00695761482616326 + 0.0285139952234768i;0.0259577390330907 - 0.0111561746553496i;-0.0234651355058812 - 0.0134231077726961i;0.00322045184215660 + 0.0255023924126257i;0.0174191236356797 - 0.0169192832056161i;-0.0225092324922808 - 0.00350639052627144i;0.0716553554322047 + 0.134175266247107i];

w4=chebwin(eleNum,9.5).*(exp(-1i*k*d4*cosd(thetaM)).');
%% array distributions
figure(1);
subplot(4,1,1),plot(d1,zeros(1,length(d1)),'kx','LineWidth',2);
title('16-element non-uniformly spaced array');
axis([0,d1(length(d1)),-0.5,0.5]);
subplot(4,1,2),plot(d2,zeros(1,length(d2)),'bo','LineWidth',2);
title('16-element \lambda/2 uniformly spaced array');
axis([0,d1(length(d1)),-0.5,0.5]);
subplot(4,1,3),plot(d3,zeros(1,length(d3)),'r^','LineWidth',2);
title('16-element uniformly spaced sparse array');
axis([0,d1(length(d1)),-0.5,0.5]);
subplot(4,1,4),plot(d4,zeros(1,length(d4)),'gs','LineWidth',2);
title('26-element \lambda/2 uniformly spaced array');
axis([0,d1(length(d1)),-0.5,0.5]);
xlabel('Wavelength');

%% Array factor
A1=zeros(length(theta),elementNumber);
for nn=1:length(d1)
    A1(:,nn)=exp(1i*k*d1(nn)*cosd(theta));
end

A2=zeros(length(theta),elementNumber);
for nn=1:length(d2)
    A2(:,nn)=exp(1i*k*d2(nn)*cosd(theta));
end

A3=zeros(length(theta),elementNumber);
for nn=1:length(d3)
    A3(:,nn)=exp(1i*k*d3(nn)*cosd(theta));
end

A4=zeros(length(theta),eleNum);
for nn=1:length(d4)
    A4(:,nn)=exp(1i*k*d4(nn)*cosd(theta));
end

%% Plot result
theta1=-90:angleStep:360-90;
figure(2);
plot(theta1,20*log10(abs(A1*w1))-max(20*log10(abs(A1*w1))),'k-','LineWidth',1.5);
hold on;
plot(theta1,20*log10(abs(A2*w2))-max(20*log10(abs(A2*w2))),'b--','LineWidth',1.5);
plot(theta1,20*log10(abs(A3*w3))-max(20*log10(abs(A3*w3))),'r:','LineWidth',1.5);
plot(theta1,20*log10(abs(A4*w4))-max(20*log10(abs(A4*w4))),'g-.','LineWidth',1.5);
hold off;
axis([-90,90,-30,0]);
xlabel('Zenith angle (Degree)');
ylabel('Normalized directivity (dBi)');
legend('16-element non-uniformly spaced','16-element \lambda/2 uniformly spaced','16-element uniformly spaced sparse array','26-element \lambda/2 uniformly spaced');
