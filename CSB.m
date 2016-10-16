clear all
close all
clc

%% INPUT
fluid='r134a';
Tcella=-10; %�C
Tamb=20;
Qev=0.2; %kW
DT_sur=0;
DT_sub=0;
DT_min_ev=10;
DT_min_covett=0:5:50;
sono bellissimo
%funzione per il rendimento del compressore f(x)=a*exp(-b*x)+c
a=-2.648;
b=1.553;
c=0.6085;
%% PUNTO 1
Tev=Tcella-DT_min_ev-DT_sur;
Pev=refpropm('p','t',Tev+273.15,'q',1,fluid)/100;
P1=Pev;
if DT_sur==0
    h1=refpropm('h','p',Pev*100,'q',1,fluid)/1000; %kJ |moltiplica per 100 perch� lo vuole in kPa
    s1=refpropm('s','p',Pev*100,'q',1,fluid)/1000; %kJ/K
else
    h1=refpropm('H','T',T1+273.15,'P',Pev*100,fluid)/1000; %kJ |moltiplica per 100 perch� lo vuole in kPa
    s1=refpropm('S','T',T1+273.15,'P',Pev*100,fluid)/1000; %kJ/K
end
T1=Tev;

%% PUNTO 3
for i=1:length(DT_min_covett)
    DT_min_co=DT_min_covett(i);
    
    T3(i)=Tamb+DT_min_co;
    Tco=T3(i)+DT_sub;
    
    Pco=refpropm('p','T',Tco+273.15,'q',0,fluid)/100;   %bar
    p3=Pco;
    if DT_sub==0
        h3=refpropm('h','p',Pco*100,'q',0,fluid)/1000;
        s3(i)=refpropm('s','p',Pco*100,'q',0,fluid)/1000;
    else
        h3=refpropm('H','T',T3(i)+273.15,'P',Pco*100,fluid)/1000;
        s3(i)=refpropm('s','T',T3(i)+273.15,'P',Pco*100,fluid)/1000;
    end
      
    s3=s3(i);
    %calcolo il rapporto di compressione
    beta=Pco/Pev;
    
    %calcolo il rendimento del compressore tramite la funzione bella
    eta=a*exp(-b*beta)+c;
    
%% PUNTO 2S
P2=Pco;
h2s=refpropm('h','p',Pco*100,'s',s1*1000,fluid)/1000; % kJ/kg
s2s=refpropm('s','p',Pco*100,'h',h2s*1000,fluid)/1000;
t2s=refpropm('t','p',Pco*100,'s',s2s*1000,fluid)-273.15;

slon=s3:0.001:s2s;
for j=1:length(slon)
    tmatr(j)=refpropm('t','p',Pco*100,'s',slon(j)*1000,fluid)-273.15;
    t2pincopallino=refpropm('t','p',Pco*100,'s',s2s*1000,fluid)-273.15;
     t3pincopallino=refpropm('t','p',Pco*100,'s',s3*1000,fluid)-273.15;
end
plot(slon,tmatr)

%% PUNTO 2
h2=h1+(h2s-h1)/eta;
t2(i)=refpropm('t','h',h2*1000,'p',Pco*100,fluid)-273.15;
s2(i)=refpropm('s','t',t2(i)+273.15,'p',Pco*100,fluid)/1000;

%% PUNTO 4
h4=h3;
P4=Pev;
t4=Tev;
s4=refpropm('s','p',P4*100,'h',h4*1000,fluid)/1000;

%% CALCOLI
m=Qev/(h1-h4); %kg/s
L=m*(h2-h1); %kj/s lavoro reale
COP(i)=Qev/L;

Exdiscp(i)=m*(s2(i)-s1)*(Tamb+273.15);
Exdisva(i)=m*(s4-s3)*(Tamb+273.15);
Exdisco(i)=m*(s3-s2(i))*(Tamb+273.15)+(Qev+L); %Qev+L: bilancio di prima specie su tutto il sistema. Sarebbe h2-h3.
Exdisev(i)=m*(s1-s4)*(Tamb+273.15)-Qev/(Tcella+273.15)*(Tamb+273.15);

Exdistot(i)=Exdiscp(i)+Exdisva(i)+Exdisco(i)+Exdisev(i);

figure(4)
plot(slon,tmatr)

end

%% CICLO T-S
tcr=refpropm('t','c',0,' ',0,fluid);
tvett=Tev+273.15-10:tcr;

for i=1:length(tvett-1)
   sls(i)=refpropm('s','t',tvett(i),'q',0,fluid);
   svs(i)=refpropm('s','t',tvett(i),'q',1,fluid);
end
sls(length(tvett))=refpropm('s','c',0,'',0,fluid);
svs(length(tvett))=sls(length(tvett));
%% PLOT
figure(1)
plot(DT_min_covett,Exdisco,'ro-')
hold on
plot(DT_min_covett,Exdiscp,'yo-')
hold on
plot(DT_min_covett,Exdisva,'ko-')
hold on
plot(DT_min_covett,Exdistot,'go-')
legend({'Exergia distrutta al condensatore','Exergia distrutta al compres.','Exergia distrutta alla valvola','Exergia distrutta totale'});

figure(2)
plot(DT_min_covett,COP,'bo-')
legend({'COP'});

figure(3)
plot(sls,tvett,svs,tvett,s2*1000,t2+273.15,'x',s3*1000,T3+273.15,'o',slon(j)*1000,tmatr(j)+273.15);
legend({'Curva inferiore','Curva superiore','Punto 2 al variare della Temperatura di Condensazione','Punto 3 al variare della Temperatura di Condensazione'})

a=polyfit(DT_min_covett,Exdisco,1);
b=polyfit(DT_min_covett,Exdistot,1);

CSB1=a(1)/b(1)


