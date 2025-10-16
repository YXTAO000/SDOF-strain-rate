clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code provides an example that can be run directly
% Cases of steel and RC specimens are given separately
% Replace the Initial parameters with your target specimen

% For more information, please see: 
% Tao Y, Huang Y, Huang X-H, Zhang B, Liu W. 
% Dynamic assessment of RC and steel structures under progressive collapse with consideration of strain rate. 
% Struct. 2025;82:110438. https://doi.org/10.1016/j.istruc.2025.110438. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Initial parameters for steel specimen
u0 = 0;
v0 = 0;
u1 = 25.2; %mm
u2 = 265.3;
u3 = 349.1;
R1 = 53090; %N
R2 = 90970;
R3 = 161100;

tr=32/1000;

me=2.5; %t
P0=35900; %N
zn=0.02;
L=2350; %mm
h=146; %mm
dt=tr/100;
Dy=40.4;
ny=5;
Dy2=6844;
ny2=3.91;

spm = 1; % Steel: 1; RC: 2


% %%%Initial parameters for RC specimen
% u0 = 0;
% v0 = 0;
% u1 = 40;
% u2 = 300;
% u3 = 610;
% R1=52100;
% R2 = 39800;
% R3 = 126200;
% %%%荷载时间
% tr=33/1000; %637ms 6�?
% %%%其余参数
% me= 415.5/1000;
% P0=  12.7*1000;
% zn=0.02;
% L=2600; 
% dt=tr/100;
% Dy=40.4;
% ny=5;
% Dy2=6844;
% ny2=3.91;
% h=250;
% D=2.57443E-05;
% 
% b=150;
% as=30; 
% At=2*16*16*3.14/4;
% Ab=2*12*12*3.14/4;
% T1=At*450/1000;
% T2=Ab*450/1000;
% C1=Ab*450/1000;
% C2=At*450/1000;
% E=200;
% mcd=40;
% fc=0.0269;
% b1=min(0.85-(fc*1000-27.58)/6.895*0.05,0.85);
% Ec=4.7*(fc*1000)^0.5;
% S=16;
% nc=(T1-T2-C1+C2)/(1.7*b1*fc*b);
% kc=((1/h/b/Ec+1/S/L)*(0.85*fc*b1*b*(h/2-mcd/4-nc)+C2-T2))/(mcd+0.85*fc*b1*0.46*L^2*b*(1/h/b/Ec+1/S/L));
% ybc=kc*mcd;
% 
% spm = 2;

%%%%%%%%%%%%k1%%%%%%%%%%%%
ut=0;
t=0;
tei=0;
count = 0;
x = zeros([],1) ;
y = zeros([],1) ;
uta=0;
R10=R1;
while ut < u1
    if spm >= 2
    nc=(At-Ab)*(h-2*as)*ut/(0.2125*fc*b1*b*L^2/E+2*At*ut+2*Ab*ut);
    c1=h/2-L^2*kc/2+nc;
    c2=h/2-L^2*kc/2-nc;
    ybsc1=8*ut/L^2*(c1-as);
    ybsc2=8*ut/L^2*(c2-as);
    ybst1=8*ut/L^2*(h-c1-as);
    ybst2=8*ut/L^2*(h-c2-as);
    yblst1=8*(uta-ut)/L^2*(h-c1-as)/dt;
    yblst2=8*(uta-ut)/L^2*(h-c2-as)/dt;
    yblsc1=8*(uta-ut)/L^2*(c1-as)/dt;
    yblsc2=8*(uta-ut)/L^2*(c2-as)/dt;
    yblc=kc*(uta-ut)/dt;

    DIFc=max((yblc/30e-6).^0.014,1);
    DIFst1=max(1+(yblst1/Dy)^(1/ny),1);
    DIFst2=max(1+(yblst2/Dy)^(1/ny),1);
    DIFsc1=max(1+(yblsc1/Dy)^(1/ny),1);
    DIFsc2=max(1+(yblsc2/Dy)^(1/ny),1);
    Rcaa1=1/(0.46*L)*(0.85*(fc*DIFc)*b1*b*h*(h/2*(1-b1/2)+mcd/4*(b1-3)+0.46*(L^2)*ybc*(b1-1)/mcd+mcd^2/8/h*(2-b1/2)+0.46*L^2*ybc*(1-b1/2)/h-b1*0.46^2*L^4*ybc^2/h/mcd^2))*1000;
    Rcaa2=1/(0.46*L)*(((C1*DIFsc1)+(C2*DIFsc2))*(h/2-as-mcd/2)+(T1*DIFst1+T2*DIFst2)*(h/2-as+mcd/2))*1000;
    Rcaa3=1/(0.46*L)*(-(T1*DIFst1-T2*DIFst2-C1*DIFsc1+C2*DIFsc2)^2/3.4/(fc*DIFc)/b)*1000;
    R1=Rcaa1+Rcaa2+Rcaa3;
    else
    ybl=3*h*(uta-ut)/dt/L^2;
    DIF=1+(ybl/Dy)^(1/ny);
    R1=R10*DIF;
    end


    a=P0/(me*tr);
    w1=(R1/(me*u1))^0.5;
    wd1=w1*(1-zn^2)^0.5;
    count = count+1 ;
    pp1=w1^2*u0;
    Ati=2*zn*a/w1^3+u0-pp1/w1^2;
    Bti=(2*zn^2-1)*a/(w1^2*wd1)+(zn*w1*u0+v0)/wd1-zn*pp1/w1/wd1;
    Cti=-2*zn*a/w1^3;
    Dti=(1-2*zn^2)*a/w1^2/wd1;
    Eti=a/w1^2;
    ut=exp(-zn*w1*(t-tei))*(Ati*cos(wd1*(t-tei))+Bti*sin(wd1*(t-tei)))+exp(-zn*w1*(t-tr))*(Cti*cos(wd1*(t-tr))+Dti*sin(wd1*(t-tr)))*hs(t,tr)+Eti*((t-tei)-(t-tr)*hs(t,tr))+Cti*(1-hs(t,tr))+pp1/w1^2;
    Atvi=(-zn*w1*Ati+wd1*Bti);
    Btvi=(-wd1*Ati-zn*w1*Bti);
    Ctvi=(-zn*w1*Cti+wd1*Dti);
    Dtvi=(-wd1*Cti-zn*w1*Dti);
    vt=exp(-zn*w1*(t-tei))*(Atvi*cos(wd1*(t-tei))+Btvi*sin(wd1*(t-tei)))+exp(-zn*w1*(t-tr))*(Ctvi*cos(wd1*(t-tr))+Dtvi*sin(wd1*(t-tr)))*hs(t,tr)+Eti*(1-hs(t,tr));
    if abs(uta-ut) > 0.0001
    uta=ut+dt*ut; 
    end
   if vt < -0.001
      disp('safe');
      disp('end-k1');
      break
   end    

x(count)=t*1000;
y(count)=real(-ut);
% y2(count)=real(vt);
vdata(count)=real(vt);
DIFdata(count)=real(DIF);  % DIFc
% ybldata(count)=real(ybl);
t=t+dt;
  
end
plot(x,y);
title('k1');
%%%initial conditions for the next stage%%%
t1=t-dt;
v1=vt;

if vt > 0
%%%%%%%%%%%%k2%%%%%%%%%%%%
k2=(R2-R10)/(u2-u1);
w2=(abs(k2)/me)^0.5;
o1=R10/me;

up=0;
t=t1;
countp = 0;
xp = zeros([],1) ;
yp = zeros([],1) ;
uta=u1;
R20=R2;
while up < u2
    if k2 < 0
         k2=(R20-R10)/(u2-u1);
    else
    ybl=3*h*(uta-up)/dt/L^2;
    DIF=1+(ybl/Dy)^(1/ny);
    R2=R20*DIF;
    k2=(R2-R1)/(u2-u1);
    end
    w2=(abs(k2)/me)^0.5;
countp = countp+1 ; 
tpi=t1;
%%Case 1 k2 > 0 (plastic hardening)
wd2=w2*(1-zn^2)^0.5;
if k2 > 0
  if t1<tr
    pp2= a*t1+w2^2*u1-R2/me;  
    Ai=2*zn*a/w2^3+u1-pp2/w2^2;
    Bi=(2*zn^2-1)*a/(w2^2*wd2)+(zn*w2*u1+v1)/wd2-zn*pp2/w2/wd2;
    Ci=-2*zn*a/w2^3;
    Di=(1-2*zn^2)*a/w2^2/wd2;
    Ei=a/w2^2;
    up=exp(-zn*w2*(t-tpi))*(Ai*cos(wd2*(t-tpi))+Bi*sin(wd2*(t-tpi)))+exp(-zn*w2*(t-tr))*(Ci*cos(wd2*(t-tr))+Di*sin(wd2*(t-tr)))*hs(t,tr)+Ei*((t-tpi)-(t-tr)*hs(t,tr))+Ci*(1-hs(t,tr))+pp2/w2^2;
    Avi=(-zn*w2*Ai+wd2*Bi);
    Bvi=(-wd2*Ai-zn*w2*Bi);
    Cvi=(-zn*w2*Ci+wd2*Di);
    Dvi=(-wd2*Ci-zn*w2*Di);
    vp=exp(-zn*w2*(t-tpi))*(Avi*cos(wd2*(t-tpi))+Bvi*sin(wd2*(t-tpi)))+exp(-zn*w2*(t-tr))*(Cvi*cos(wd2*(t-tr))+Dvi*sin(wd2*(t-tr)))*hs(t,tr)+Ei*(1-hs(t,tr));
    if abs(uta-up) > 0.0001
    uta=up+dt*up;
    end
    if vp < -0.001
      disp('safe');
      disp('end-Case1-t1<tr');
      break
   end
  else
    yy2=a*tr+w2^2*u1-R10/me;
    Fi=u1-yy2/w2^2;
    Gi=(zn*w2*u1+v1)/wd2-zn*yy2/w2/wd2;
    Fvi=(-zn*w2*Fi+wd2*Gi);
    Gvi=(-wd2*Fi-zn*w2*Gi);
    up=exp(-zn*w2*(t-tpi))*(Fi*cos(wd2*(t-tpi))+Gi*sin(wd2*(t-tpi)))+yy2/w2^2;
    vp=exp(-zn*w2*(t-tpi))*(Fvi*cos(wd2*(t-tpi))+Gvi*sin(wd2*(t-tpi)));
    if abs(uta-up) > 0.0001
    uta=up+dt*up;  %0.01
    end
    if vp < -0.001
      disp('safe');
      disp('end-Case1-t1>=tr');
      break
   end
  end
end

%%Case 2 k2 = 0 (perfect plasticity)
if k2 == 0
  if t1<tr  
    up=a/6*(t-t1)^3+(a*t1-o1)/2*(t-t1)^2+v1*(t-t1)+u1-a/6*(t-tr)^3*hs(t,tr);
    vp=a/2*(t-t1)^2+(a*t1-o1)*(t-t1)+v1-a/2*(t-tr)^2*hs(t,tr);
    if abs(uta-up) > 0.0001
    uta=up+dt*up; 
    end
    if vp < -0.001
      disp('safe');
      disp('end-Case2-t1<tr');
      break
   end
  else
    up=(a*tr-o1)*(t-t1)^2/2+v1*(t-t1)+u1;
    vp=(a*tr-o1)*(t-t1)+v1;
    if abs(uta-up) > 0.0001
    uta=up+dt*up;  
    end
    if vp < -0.001 
      disp('safe');
      disp('end-Case2-t1>=tr');
      break
   end
  end 
end

%%Case 3 k2 < 0 (plastic softening)
wb2=w2*(1+zn^2)^0.5;
if k2 < 0
  if t1<tr
    oo2=a*t1-w2^2*u1-v1;  
    Asi=2*zn*a/w2^3+u1+oo2/w2^2;
    Bsi=(2*zn^2+1)*a/(w2^2*wb2)+(zn*w2*u1+v1)/wb2+zn*oo2/w2/wb2;
    Csi=-2*zn*a/w2^3;
    Dsi=-(1+2*zn^2)*a/w2^2/wb2;
    Esi=-a/w2^2;
    Asvi=(-zn*w2*Asi+wd2*Bsi);
    Bsvi=(wd2*Asi-zn*w2*Bsi);
    Csvi=(-zn*w2*Csi+wd2*Dsi);
    Dsvi=(wd2*Csi-zn*w2*Dsi);
    up=exp(-zn*w2*(t-tpi))*(Asi*cosh(wb2*(t-tpi))+Bsi*sinh(wb2*(t-tpi)))+exp(-zn*w2*(t-tr))*(Csi*cosh(wb2*(t-tr))+Dsi*sinh(wb2*(t-tr)))*hs(t,tr)+Esi*((t-tpi)-(t-tr)*hs(t,tr))+Csi*(1-hs(t,tr))-oo2/w2^2;
    vp=exp(-zn*w2*(t-tpi))*(Asvi*cosh(wb2*(t-tpi))+Bsvi*sinh(wb2*(t-tpi)))+exp(-zn*w2*(t-tr))*(Csvi*cosh(wb2*(t-tr))+Dsvi*sinh(wb2*(t-tr)))*hs(t,tr)+Esi*(1-hs(t,tr));
    if abs(uta-up) > 0.0001
    uta=up+dt*up;
    end
    if vp < -0.001
      disp('safe');
      disp('end-Case3-t1<tr');
      break
   end
  else
    nn2=a*tr-w2^2*u1-R10/me;
    Fsi=u1+nn2/w2^2;
    Gsi=(zn*w2*u1+v1)/wb2+zn*nn2/w2/wb2;
    Fsvi=(-zn*w2*Fsi+wb2*Gsi);
    Gsvi=(wb2*Fsi-zn*w2*Gsi);
    up=exp(-zn*w2*(t-tpi))*(Fsi*cosh(wb2*(t-tpi))+Gsi*sinh(wb2*(t-tpi)))-nn2/w2^2;
    vp=exp(-zn*w2*(t-tpi))*(Fsvi*cosh(wb2*(t-tpi))+Gsvi*sinh(wb2*(t-tpi)));
    if abs(uta-up) > 0.0001
    uta=up+dt*up;
    end
    if vp < -0.001
      disp('safe');
      disp('end-Case3-t1>=tr');
      break
   end
  end
end
xp(countp)=t*1000;
yp(countp)=real(-up); 
% vdata(count+countp)=real(vp);
DIFdata(count+countp)=real(DIF);
% ybldata(count+countp)=real(ybl);
t=t+dt; 
end
plot(xp,yp);
title('k2');


xtp=[x xp];
ytp=[y yp];
plot(xtp,ytp);
title('k1-k2');

%%%initial conditions for the next stage%%%
t2=t-dt;
v2=vp;

end

if vt > 0 && vp> 0

%%%%%%%%%%%%k3%%%%%%%%%%%%
o2=R20/me;

uca=0;
t=t2;
tca=t2;
countca = 0;
xca = zeros([],1) ;
yca = zeros([],1) ;
uta=u2;
R30=R3;
while uca < u3
    ybl=((1+(uta/L)^2)^0.5-(1+(uca/L)^2)^0.5)/dt;
    DIF=1+(ybl/Dy2)^(1/ny2);
    R3=R30*DIF;
    k3=(R3-R2)/(u3-u2);
    w3=(k3/me)^0.5;
pp3= a*t2+w3^2*u2-R3/me;
wd3=w3*(1-zn^2)^0.5;
countca = countca+1;         
if t2<tr
    Ai=2*zn*a/w3^3+u2-pp3/w3^2;
    Bi=(2*zn^2-1)*a/(w3^2*wd3)+(zn*w3*u2+v2)/wd3-zn*pp3/w3/wd3;
    Ci=-2*zn*a/w3^3;
    Di=(1-2*zn^2)*a/w3^2/wd3;
    Ei=a/w3^2;
    uca=exp(-zn*w3*(t-tca))*(Ai*cos(wd3*(t-tca))+Bi*sin(wd3*(t-tca)))+exp(-zn*w3*(t-tr))*(Ci*cos(wd3*(t-tr))+Di*sin(wd3*(t-tr)))*hs(t,tr)+Ei*((t-tca)-(t-tr)*hs(t,tr))+Ci*(1-hs(t,tr))+pp3/w3^2;
    Avi=(-zn*w3*Ai+wd3*Bi);
    Bvi=(-wd3*Ai-zn*w3*Bi);
    Cvi=(-zn*w3*Ci+wd3*Di);
    Dvi=(-wd3*Ci-zn*w3*Di);
    vca=exp(-zn*w3*(t-tca))*(Avi*cos(wd3*(t-tca))+Bvi*sin(wd3*(t-tca)))+exp(-zn*w3*(t-tr))*(Cvi*cos(wd3*(t-tr))+Dvi*sin(wd3*(t-tr)))*hs(t,tr)+Ei*(1-hs(t,tr));
    if abs(uta-uca) > 0.0001
    uta=uca+dt*uca;
    end
  if vca < -0.001 
      disp('safe');
      disp('end-CA-t2<tr');
      break
   end
else
    yy3=a*tr+w3^2*u2-R2/me;
    Fi=u2-yy3/w3^2;
    Gi=(zn*w3*u2+v2)/wd3-zn*yy3/w3/wd3;
    Fvi=(-zn*w3*Fi+wd3*Gi);
    Gvi=(-wd3*Fi-zn*w3*Gi);
    uca=exp(-zn*w3*(t-tca))*(Fi*cos(wd3*(t-tca))+Gi*sin(wd3*(t-tca)))+yy3/w3^2;
    vca=exp(-zn*w3*(t-tca))*(Fvi*cos(wd3*(t-tca))+Gvi*sin(wd3*(t-tca)));
    if abs(uta-uca) > 0.0001
    uta=uca+dt*uca; 
    end
  if vca < -0.001 
      disp('safe');
      disp('end-CA-t2>=tr');
      break
   end
end
xca(countca)=t*1000;
yca(countca)=real(-uca);
% vdata(count+countp+countca)=real(vca);
DIFdata(count+countp+countca)=real(DIF);
% ybldata(count+countp+countca)=real(ybl);
t=t+dt; 
end
plot(xca,yca);
title('k3');

xtpca=[xtp xca];
ytpca=[ytp yca];
plot(xtpca,ytpca);
title('k1-k3');

if uca > u3
   disp('collapse') 
end

end