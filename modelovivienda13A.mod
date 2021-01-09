% DSGE CON POLITICA FISCAL, DEMANDA DE DINERO, MERCADO FINANCIERO, MERCADO DE VIVIENDA Y CRÉDITO HIPOTECARIO
% Colombia 2015
% Modelo base con oferta y demanda de exportaciones.
% Precios flexibles.
% Incluye Inversión extranjera directa.
% Exportaciones extractivas (petróleo), que se reparten entre empresas y gobierno, mediante el "state take".
% Competencia monopolística. La ganancia es transferida a los inversionistas para financiar la inversión.
% Demanda y oferta de dinero.
% Separación de las decisiones de ahorrar e invertir.
% El intermediario financiero maximiza su ganancia, dada una función de pérdida, costos a depósitos y crédito del banco central.
% El banco central determina la oferta monetaria.
% Los hogares obtienen su utilidad del consumo, el ocio, la propiedad de vivienda y las tenencias de dinero.
% Hay un sector productivo, que distribuye si produccion entre bienes y servicios y construccion

close all;
external_function(name=logncdf, nargs=3);

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var 
c
n
w
ir
bf
er
q
i 
ifdi
k
kfdi
kv
r
rfdi
fi
df
py
y
g
d 
m
s
x
pm
pd
px
xx
xt
pwx
pib
ri
tauh
tauk
f 
fk
z
iva
aran
dm
pi 
vr
ppet
qpet
lamda
ga
po
v1
v2
iri
gb
rb
irb
tc
cb
ira
walras
scc
sck
pw
iv
rv
pv
yv
yt
pt
taukv
kapa
ltv
irh
ch
tauv
rrb
psiv
psim
psin
u
pomeg
omeg
rs
sdomeg
chi
ct
%e24
;
varexo e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 e11 e12 e13 e14 e15 e16 e17 e18 e19 e20 e21 e22 e23 
e24 
e25 e26;

parameters 
psin0
beta
theta
tauh0
delta
gama
tauk0
f0
fk0
alfa1
alfa2
alfa3
z0
b
sig
omega
iva0
g0
x0
px0
pwx0
omegae
be
sige
omegad
bd
sigd
xt0
pw0
pwm
ma
rho0
rm
a
meta
aran0
tauhk
taukk
ivak
arank
psim0
vr0
ifdi0
ir0
gk
ppet0
qpet0
share
epsilon
rentas
rcc
oabc
irb0
rrb0
alfaf
kf
alfah
kh
iri0
sigmatc
deltatc
gamac
gaman
gamam
gamav
pib0
pi0
alfay
alfapi
kbb
vbb
rr 
spread
spread2
deltav
kapa0
tauv0
psiv0
bv
sigv
omegav
taukv0
sltv
ppch
ltv0
chi0
act
mu
sdomeg0
rs0
lamdi
lamdv

rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8, rho9, rho10 rho11 rho12 rho13 rho14 rho15 rho16 rho17 rho18 rho19 rho20  rho21 rho22 rho23 rho24 rho25 rho26 cor1;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

psin0	=	 0.839925651 	;
beta	=	0.966660992	;
theta	=	 1.011289 	;
tauh0	=	 0.13921 	;
delta	=	0.942928263	;
gama	=	1.561073496	;
tauk0	=	 0.133974 	;
f0	=	0.760360883	;
fk0	=	0.80712262	;
alfa1	=	0.276098782	;
alfa2	=	0.023635479	;
alfa3	=	0.030266205	;
z0	=	 26.59720303 	;
b	=	2.064665685	;
sig	=	1.5	;
omega	=	0.330468306	;
iva0	=	0.123953184	;
g0	=	3.946133413	;
x0	=	1.88796393	;
px0	=	1	;
pwx0	=	1	;
omegae	=	0.226246806	;
be	=	1.451223777	;
sige	=	3	;
omegad	=	0.712316809	;
bd	=	2.413222496	;
sigd	=	-2	;
xt0	=	77.40652112	;
pw0	=	1	;
pwm	=	0.972498643	;
ma	=	0.052631579	;
rho0	=	0.5	;
rm	=	 0.021307 	;
a	=	0.061189212	;
meta	=	0.408404801	;
aran0	=	0.02827907	;
tauhk	=	1	;
taukk	=	2	;
ivak	=	2	;
arank	=	2	;
psim0	=	 0.07377068 	;
vr0	=	0.036093325	;
ifdi0	=	0.247142369	;
ir0	=	 1.06552350 	;
gk	=	2	;
ppet0	=	1	;
qpet0	=	0.801901957	;
share	=	0.292740852	;
epsilon	=	20	;
rentas	=	1.832450981	;
rcc	=	0.476790463	;
oabc	=	0.037022972	;
irb0	=	 1.0570588 	;
rrb0	=	0.076836026	;
alfaf	=	1.2	;
kf	=	0.017660472	;
alfah	=	1.2	;
kh	=	0.040510778	;
iri0	=	 1.10215157 	;
sigmatc	=	 10.00000000 	;
deltatc	=	 1.44595993 	;
gamac	=	 0.57000000 	;
gaman	=	 2.50000000 	;
gamam	=	 2.00000000 	;
gamav	=	 5.00000000 	;
pib0	=	17.18740916	;
pi0	=	 1.0300000 	;
alfay	=	 0.5000000 	;
alfapi	=	 1.5000000 	;
kbb	=	-0.08470711	;
vbb	=	-0.272727273	;
rr 	=	 0.0270588 	;
spread	=	0.036851467	;
spread2	=	0.008464672	;
deltav	=	0.909317602	;
kapa0	=	1	;
tauv0	=	 1.0000000 	;
psiv0	=	404.907184	;
bv	=	 2.950827 	;
sigv	=	-2	;
omegav	=	0.796372447	;
taukv0	=	1	;
sltv	=	0.500000000	;
ppch	=	0.403621206	;
ltv0	=	0.700000000	;
chi0	=	0.001725	;
act	=	0.957405739	;
mu	=	0	;
sdomeg0	=	0.5	;
rs0	=	0.115	;
lamdi	=	0.5	;
lamdv	=	0.75	;

psin0	=	 0.839925651 	;
beta	=	0.966660992	;
theta	=	 1.011289 	;
tauh0	=	 0.13921 	;
delta	=	0.942928263	;
gama	=	1.561073496	;
tauk0	=	 0.133974 	;
f0	=	0.760360883	;
fk0	=	0.80712262	;
alfa1	=	0.276098782	;
alfa2	=	0.023635479	;
alfa3	=	0.030266205	;
z0	=	 26.59720303 	;
b	=	2.064665685	;
sig	=	1.5	;
omega	=	0.330468306	;
iva0	=	0.123953184	;
g0	=	3.946133413	;
x0	=	1.88796393	;
px0	=	1	;
pwx0	=	1	;
omegae	=	0.226246806	;
be	=	1.451223777	;
sige	=	3	;
omegad	=	0.712316809	;
bd	=	2.413222496	;
sigd	=	-2	;
xt0	=	77.40652112	;
pw0	=	1	;
pwm	=	0.972498643	;
ma	=	0.052631579	;
rho0	=	0.5	;
rm	=	 0.021307 	;
a	=	0.061189212	;
meta	=	0.408404801	;
aran0	=	0.02827907	;
tauhk	=	1	;
taukk	=	2	;
ivak	=	2	;
arank	=	2	;
psim0	=	 0.07377068 	;
vr0	=	0.036093325	;
ifdi0	=	0.247142369	;
ir0	=	 1.06552350 	;
gk	=	2	;
ppet0	=	1	;
qpet0	=	0.801901957	;
share	=	0.292740852	;
epsilon	=	20	;
rentas	=	1.832450981	;
rcc	=	0.476790463	;
oabc	=	0.037022972	;
irb0	=	 1.0570588 	;
rrb0	=	0.076836026	;
alfaf	=	1.2	;
kf	=	0.017660472	;
alfah	=	1.2	;
kh	=	0.040510778	;
iri0	=	 1.10215157 	;
sigmatc	=	 10.00000000 	;
deltatc	=	 1.44595993 	;
gamac	=	 0.57000000 	;
gaman	=	 2.50000000 	;
gamam	=	 2.00000000 	;
gamav	=	 5.00000000 	;
pib0	=	17.18740916	;
pi0	=	 1.0300000 	;
alfay	=	 0.5000000 	;
alfapi	=	 1.5000000 	;
kbb	=	-0.08470711	;
vbb	=	-0.272727273	;
rr 	=	 0.0270588 	;
spread	=	0.036851467	;
spread2	=	0.008464672	;
deltav	=	0.909317602	;
kapa0	=	1	;
tauv0	=	 1.0000000 	;
psiv0	=	404.907184	;
bv	=	 2.950827 	;
sigv	=	-2	;
omegav	=	0.796372447	;
taukv0	=	1	;
sltv	=	0.500000000	;
ppch	=	0.403621206	;
ltv0	=	0.700000000	;
chi0	=	0.001725	;
act	=	0.957405739	;
mu	=	0	;
sdomeg0	=	0.5	;
rs0	=	0.115	;
lamdi	=	0.5	;
lamdv	=	0.75	;

    
rho1    = 0.95;
rho2    = 0.95;
rho3    = 0.95;
rho4    = 0.95;
rho5    = 0.95;
rho6    = 0.95;
rho7    = 0.95;
rho8    = 0.95;
rho9    = 0.95;
rho10   = 0.95;
rho11   = 0.95;
rho12   = 0.95;
rho13   = 0.95;
rho14   = 0.95;
rho15   = 0.95;
rho16   = 0.95;
rho17   = 0.95;
rho18   = 0.95;
rho19   = 0.95;
rho20   = 0.95;
rho21   = 0.95;
rho22   = 0.95;
rho23   = 0.95;
rho24   = 0.95;
rho25   = 0.95;
rho26   = 0.95;
cor1    = 0.95;
% sigma   = 0.2222;
sigma=0.01;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 

% Modelación de los hogares (1 a 7)
  (exp(c)^gamac/(1-exp(n))^gaman)=(1-exp(tauh))*exp(w)*exp(kapa)/exp(psin);
  (exp(c(+1))/exp(c))^gamac=beta*exp(ir(+1))/exp(pi(+1));
  exp(dm)=(exp(c)^(gamac/gamam))*(exp(kapa)*(exp(ir(+1))-1)/(exp(ir(+1))*exp(psim)))^(-1/gamam);
  exp(kv)=((1/exp(psiv))*(exp(kapa)*(exp(c))^(-gamac)*(exp(taukv)*exp(pv)*(1-exp(ltv)*sltv))-beta*exp(kapa(+1))*(exp(c(+1)))^(-gamac)*(exp(tauv(+1))*exp(rv(+1))
             +exp(taukv(+1))*exp(pv(+1))*deltav*(1-exp(ltv(+1))*sltv))))^(-1/gamav);
  exp(bf)-exp(bf(-1))*exp(ir)/(theta*exp(pi))=exp(ch)-exp(ch(-1))/(theta*exp(pi))+(1-exp(tauh))*(exp(w)*exp(n)+exp(ga))
             +(exp(tauv))*exp(rv)*exp(kv(-1))/theta+rentas+exp(gb)+exp(er)*(exp(f)+exp(fk))+exp(dm(-1))/(theta*exp(pi))
             -exp(c)-exp(dm)-exp(taukv)*exp(pv)*exp(iv)+exp(ltv)*sltv*exp(taukv)*exp(pv)*exp(iv)-ppch*exp(ch(-1))/(theta*exp(pi))
             -(exp(irh)-1)*exp(ch(-1))/(theta*exp(pi));   
  exp(ch)=exp(ch(-1))*exp(irh)/(theta*exp(pi))-ppch*exp(ch(-1))/(theta*exp(pi))+sltv*exp(ltv)*exp(taukv)*exp(pv)*exp(iv);
  exp(kv)=deltav*exp(kv(-1))/theta+exp(iv);

% Modelación de la inversión (8 a 12)
  exp(q)=(1/gama)*(exp(fi(+1))*exp(pi(+1))/exp(iri(+1))-1); 
  exp(q)=exp(i)/exp(k(-1));
  exp(fi)=(1-exp(tauk))*exp(r)-(gama/2)*(exp(q))^2-exp(q)+exp(fi(+1))*(delta+exp(q))*exp(pi(+1))/exp(iri(+1));
  exp(df)=exp(i)+rentas+exp(df(-1))*exp(iri)/(theta*exp(pi))-(1-exp(tauk))*(exp(r)*exp(k(-1))/theta)-(1-share)*exp(ppet)*exp(qpet)*exp(er);
  exp(k)=delta*exp(k(-1))/theta+exp(i);

% Modelación del sector bancario  (13 a 19)
  exp(ct)=exp(df)+exp(ch);
  exp(cb)/exp(bf)=(deltatc*(exp(irb(+1))/exp(ir(+1))))^(-sigmatc);  
  exp(tc(-1))*exp(ira)=exp(bf(-1))*exp(ir)+exp(cb(-1))*exp(irb);
  exp(tc)=exp(bf)+exp(cb);
  exp(rb)=exp(rrb)*(rcc*exp(dm)+exp(tc));
  rcc*exp(dm)+exp(tc)=exp(df)+exp(rb)+exp(ch);
  exp(gb)=exp(df(-1))*exp(iri)/(theta*exp(pi))+exp(ch(-1))*exp(irh)/(theta*exp(pi))+exp(rb(-1))/(theta*exp(pi))
             +exp(bf)+exp(cb)+rcc*exp(dm)-exp(bf(-1))*exp(ir)/(theta*exp(pi))-exp(cb(-1))*exp(irb)/(theta*exp(pi))
             -exp(df)-exp(ch)-exp(rb)-rcc*exp(dm(-1))/(theta*exp(pi))
             +ppch*exp(ch(-1))/(theta*exp(pi))-sltv*exp(ltv)*exp(taukv)*exp(pv)*exp(iv);

% Modelación del riesgo  (20 a 23) ajuste parcial. Para evitarlo, hacer lamdi y lamdv igual a 0
  exp(iri)=lamdi*(exp(iri(-1)))+(1-lamdi)*(-exp(rrb)/(1-exp(rrb))+exp(ira)/(1-exp(rrb))+kf*alfaf*exp(df)^(alfaf-1)+exp(chi)*exp(pomeg));
  exp(irh)=lamdv*(exp(irh(-1)))+(1-lamdv)*(-exp(rrb)/(1-exp(rrb))+exp(ira)/(1-exp(rrb))+kh*alfah*exp(ch)^(alfah-1)+exp(chi)*exp(pomeg));
%  exp(omeg)=(((exp(bf(-1))*exp(ir))+(exp(cb(-1))*exp(irb))+rcc*exp(dm(-1))-exp(rb(-1))-act-(1-exp(rs))*exp(irh)*exp(ch(-1))))/((1-exp(rs))*(exp(df(-1))*exp(iri)));  
  exp(omeg)=(exp(rs)*(exp(df(-1))*exp(iri)+exp(ch(-1))*exp(irh))-exp(rb(-1))-act+exp(bf(-1))*exp(ir)+
             exp(cb(-1))*exp(irb)+rcc*exp(dm(-1)))/(exp(df(-1))*exp(iri)+exp(ch(-1))*exp(irh));
  exp(pomeg)=logncdf(exp(omeg),mu,exp(sdomeg));

% Modelación de la inversión extranjera directa  (24).
  exp(kfdi)=delta*exp(kfdi(-1))/theta+exp(ifdi)*exp(er);

% Modelación de las empresas (25 a 34)
  alfa1*(exp(k(-1))^(alfa1-1))*(exp(kfdi(-1))^(alfa2))*(exp(kv(-1))^(alfa3))*((theta)^(1-alfa1-alfa2-alfa3))*(exp(z)*exp(n))^(1-alfa1-alfa2-alfa3)=
             exp(r)/exp(lamda);
  alfa2*(exp(kfdi(-1))^(alfa2-1))*(exp(k(-1))^(alfa1))*(exp(kv(-1))^(alfa3))*((theta)^(1-alfa1-alfa2-alfa3))*(exp(z)*exp(n))^(1-alfa1-alfa2-alfa3)=
             exp(rfdi)/exp(lamda);
  alfa3*(exp(kv(-1))^(alfa3-1))*(exp(k(-1))^(alfa1))*(exp(kfdi(-1))^(alfa2))*((theta)^(1-alfa1-alfa2-alfa3))*(exp(z)*exp(n))^(1-alfa1-alfa2-alfa3)=
             exp(rv)/exp(lamda);
  (1-alfa1-alfa2-alfa3)*exp(k(-1))^alfa1*exp(kfdi(-1))^alfa2*exp(kv(-1))^alfa3*theta^(-alfa1-alfa2-alfa3)*exp(z)^(1-alfa1-alfa2-alfa3)*exp(n)^(-alfa1-alfa2-alfa3)=
             exp(w)/exp(lamda);
  exp(y)=((exp(k(-1))/(theta))^alfa1)*((exp(kfdi(-1))/(theta))^alfa2)*((exp(kv(-1))/(theta))^alfa3)*(exp(z)*exp(n))^(1-alfa1-alfa2-alfa3);
  exp(po)=(1+ma)*(exp(v1)/exp(v2));
  exp(v1)=exp(y)*exp(lamda)+beta*rho0*exp(v1(+1));
  exp(v2)=exp(y)+beta*rho0*exp(v2(+1));
  exp(py)=(((1-rho0)*(exp(po))^(1-epsilon)+(rho0)*(exp(py(-1)))^(1-epsilon))^(1/(1-epsilon)));
  exp(ga)=exp(y)*(exp(py)-exp(lamda));

% Modelación de la frontera de produccion vivienda y otros bienes  (35 a 38)
  exp(y)=bv*(omegav*exp(yv)^((sigv-1)/sigv)+(1-omegav)*exp(yt)^((sigv-1)/sigv))^(sigv/(sigv-1));
  exp(yv)/exp(yt)=(exp(pv)*(1-omegav)/(exp(pt)*omegav))^(-sigv);
  exp(py)*exp(y)=exp(pv)*exp(yv)+exp(pt)*exp(yt);
  exp(yv)=exp(iv);

% Modelación de la demanda de importaciones (39 a 42)
  exp(c)+exp(i)+exp(ifdi)*exp(er)+exp(g)=b*(omega*exp(m)^((sig-1)/sig)+(1-omega)*exp(d)^((sig-1)/sig))^(sig/(sig-1));
  exp(m)/exp(d)=(exp(pm)*(1-omega)/(exp(pd)*omega))^(-sig);
  exp(pm)=exp(er)*pwm*(1+exp(aran));
  exp(c)+exp(i)+exp(ifdi)*exp(er)+exp(g)=(1+exp(iva))*(exp(pm)*exp(m)+exp(pd)*exp(d));

% Modelación de la demanda de exportaciones   (43 a 44)
  exp(xt)=be*(omegae*exp(x)^((sige-1)/sige)+(1-omegae)*exp(xx)^((sige-1)/sige))^(sige/(sige-1));
  exp(x)/exp(xx)=(exp(pwx)*(1-omegae)/(exp(pw)*omegae))^(-sige);

% Modelación de la oferta de exportaciones  (45 a 48)
  exp(yt)=bd*(omegad*exp(x)^((sigd-1)/sigd)+(1-omegad)*exp(d)^((sigd-1)/sigd))^(sigd/(sigd-1));
  exp(x)/exp(d)=(exp(px)*(1-omegad)/(exp(pd)*omegad))^(-sigd);
  exp(pt)*exp(yt)=exp(px)*exp(x)+exp(pd)*exp(d);
  exp(px)=exp(er)*exp(pwx);

% Restricción presupuestal del gobierno    (49 a 51)
  exp(s)*exp(er)=(((exp(er)*exp(s(-1)))*(1+exp(ri)))/(theta))-exp(tauk)*(((exp(r)*exp(k(-1)))/(theta)))-exp(tauh)*(exp(w)*exp(n)+exp(ga))
               -(exp(tauv)-1)*exp(rv)*exp(kv(-1))/theta-(exp(iva)*(exp(c)+exp(i)+exp(ifdi)*exp(er)+exp(g)))/(1+exp(iva))
               -exp(aran)*pwm*exp(er)*exp(m)-(exp(taukv)-1)*exp(pv)*exp(iv)+exp(g)-share*exp(ppet)*exp(qpet)*exp(er)-oabc;
  exp(ri)=rm+a*exp(s)*exp(er)/exp(pib);
  exp(pib)=(exp(c)+exp(g)+exp(i)+exp(iv)*exp(pv))+exp(x)+exp(ifdi)*exp(er)-exp(m)/(1+exp(aran))+exp(qpet); 

% Cuenta corriente de la balanza de pagos   (52 a 54)
  walras=exp(pwx)*exp(x)+exp(f)+exp(fk)+exp(ppet)*exp(qpet)+exp(s)+exp(ifdi)-(vr)-(((pwm)*exp(m)))-
                (((exp(s(-1)))*(1+exp(ri)))/(theta))-exp(kfdi(-1))*exp(rfdi)/(theta*exp(er));
  scc=(exp(pwx)*exp(x)-(((pwm)*exp(m))))/(exp(pib)/exp(er));
  sck=(exp(f)+exp(ppet)*exp(qpet)+exp(fk)+exp(ifdi)+exp(s)-exp(s(-1))/theta-(((exp(s(-1)))*(exp(ri)))/(theta))-
                 exp(kfdi(-1))*exp(rfdi)/(theta*exp(er)))/(exp(pib)/exp(er));

% Oferta monetaria (55)
  exp(dm)*(1-rcc)=exp(dm(-1))*(1-rcc)/(theta*exp(pi))+(vr)*exp(er)+oabc-exp(rb)+exp(rb(-1))/(theta*exp(pi))+exp(cb)-exp(cb(-1))*exp(irb)/(theta*exp(pi));

% Evaluación del bienestar (56)
  exp(u)=exp(kapa)*exp(c)^(1-gamac)/(1-gamac)+exp(psin)*(1-exp(n))^(1-gaman)/(1-gaman)+exp(psiv)*exp(kv)^(1-gamav)/(1-gamav)+exp(psim)*exp(dm)^(1-gamam)/(1-gamam)+beta*theta*exp(u(+1));

  
% Procesos exógenos  (57 a 82)
  exp(z)=(exp(z(-1))^rho1)*(z0^(1-rho1))*exp(e1);
  exp(f)=(exp(f(-1))^rho2)*(f0^(1-rho2))*exp(e2);
  exp(ifdi)=(exp(ifdi(-1))^rho3)*(ifdi0^(1-rho3))*exp(e3);
  exp(ppet)=(exp(ppet(-1))^rho4)*(ppet0^(1-rho4))/exp(e4);
  exp(qpet)=(exp(qpet(-1))^rho5)*(qpet0^(1-rho5))*exp(e5);
  exp(pw)=(exp(pw(-1))^rho6)*(pw0^(1-rho6))*exp(e6);
  exp(fk)=(exp(fk(-1))^rho7)*(fk0^(1-rho7))*exp(e7);
 %exp(g)=(exp(g(-1))^rho8)*(g0^(1-rho8))*exp(e8);
 exp(tauh)=(exp(tauh(-1))^rho9)*(tauh0^(1-rho9))*exp(e9);
  exp(iva)=(exp(iva(-1))^rho10)*(iva0^(1-rho10))*exp(e10); 
  exp(aran)=(exp(aran(-1))^rho11)*(aran0^(1-rho11))*exp(e11); 
 exp(tauk)=(exp(tauk(-1))^rho12)*(tauk0^(1-rho12))*exp(e12);
  exp(irb)=(exp(irb(-1))^rho13)*(irb0^(1-rho13))*exp(e13);
  vr =((vr(-1))^rho14)*(vr0^(1-rho14))*exp(e14);
  exp(xt)=(exp(xt(-1))^rho15)*(xt0^(1-rho15))*exp(e15);
  exp(taukv)=((exp(taukv(-1)))^rho16)*(taukv0^(1-rho16))/exp(e16);  
  exp(kapa)=((exp(kapa(-1)))^rho17)*(kapa0^(1-rho17))/exp(e17);
  exp(ltv)=((exp(ltv(-1)))^rho18)*(ltv0^(1-rho18))*exp(e18);
  exp(tauv)=((exp(tauv(-1)))^rho19)*(tauv0^(1-rho19))*exp(e19);  
  exp(rrb)=((exp(rrb(-1)))^rho20)*(rrb0^(1-rho20))*exp(e20);  
  exp(psin)=((exp(psin(-1)))^rho21)*(psin0^(1-rho21))*exp(e21);  
  exp(psim)=((exp(psim(-1)))^rho22)*(psim0^(1-rho22))*exp(e22);  
  %exp(psiv)=((exp(psiv(-1)))^rho23)*(psiv0^(1-rho23))/exp(e23);  
  exp(psiv)=((exp(psiv(-1)))^rho23)*(psiv0^(1-rho23))/exp(e23); 
  exp(rs)=(exp(rs(-1))^rho24)*(rs0^(1-rho24))*exp(e24);
  exp(sdomeg)=(exp(sdomeg(-1))^rho25)*(sdomeg0^(1-rho25))*exp(e25);
  exp(chi)=(exp(chi(-1))^rho26)*(chi0^(1-rho26))*exp(e26);

%exp(e24)=exp((e18)^(1.5));
 
% Cierre fiscal del modelo (desbloquear una y bxloquear el correspondiente proceso exógeno) (52) 
%exp(tauh)-tauh0-tauhk*(exp(s(-1))*exp(er(-1))/exp(pib(-1))-meta);
 %exp(tauk)-tauk0-taukk*(exp(s(-1))*exp(er(-1))/exp(pib(-1))-meta);
%  exp(iva(+1))-iva0-ivak*(exp(s)*exp(er)/exp(pib)-meta);
  %exp(iva)-iva0-ivak*(exp(s(-1))*exp(er(-1))/exp(pib(-1))-meta);
%  exp(aran)-aran0-arank*(exp(s(-1))*exp(er)/exp(pib(-1))-meta);
 exp(g)=g0-gk*(exp(s(-1))*exp(er(-1))/exp(pib(-1))-meta);
 
% ecuacion para cierre alternativo
%  exp(s)*exp(er)/exp(pib)=meta;

% regla de Taylor (si se activa, debe bloquearse el proceso exógeno de irb)
% exp(irb)=irb0+alfay*(exp(pib(-1))/pib0-1)+alfapi*(exp(pi(-1))/pi0-1);
% exp(irb)=irb0+alfay*(exp(pib)/pib0-1)+alfapi*(exp(pi)/pi0-1);  
% exp(irb)=irb0+alfapi*(exp(pi)/pi0-1);
% exp(irb(+1))=irb0+alfay*(exp(pib)/pib0-1)+alfapi*(exp(pi)/pi0-1);
% exp(irb)=irb0+0*(exp(pib(-1))/pib0-1)+alfapi*(exp(pi(-1))/pi0-1);
% exp(irb)=rr+exp(pi)+0.5*(exp(pi)-pi0)+0.5*(exp(pib)/pib0-1);
% exp(irb)=rr+exp(pi)+1.353*(exp(pi)-pi0)+0.064*(exp(pib)/pib0-1);

% politica cambiaria (si se activa, debe bloquearse el proceso exógeno de vr)
% vr=(sck+kbb+vbb*sck)*(exp(pib)/exp(er));

% politica fiscal contracíclica
% exp(g)=g0+pib0-exp(pib);

%Perturbación compuesta
%  exp(e1)=(exp(e2(-3))^eta)*exp(e0);


%exp(irpi)=exp(ir)/exp(pi);
%exp(irapi)=exp(ira)/exp(pi);
%exp(iripi)=exp(iri)/exp(pi);


 end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
c	=	2.459706244	;
n	=	-0.754458899	;
w	=	2.96594119	;
ir	=	0.063466224	;
bf	=	1.897751157	;
er	=	0	;
q	=	-2.694181374	;
i 	=	0.486164058	;
ifdi	=	-1.397790717	;
k	=	3.180345432	;
kfdi	=	1.296390656	;
kv	=	2.166662188	;
r	=	-1.84415562	;
rfdi	=	-2.418210605	;
fi	=	0.16802575	;
df	=	1.842207949	;
py	=	0	;
y	=	2.663253848	;
g	=	1.372736217	;
d 	=	2.448824812	;
m	=	1.389723192	;
s	=	1.948680653	;
x	=	0.635498962	;
pm	=	0	;
pd	=	0	;
px	=	0	;
xx	=	4.324378417	;
xt	=	4.349071029	;
pwx	=	0	;
pib	=	2.84417709	;
ri	=	-3.072678803	;
tauh	=	-1.971772928	;
tauk	=	-2.010112603	;
f 	=	-0.273962113	;
fk	=	-0.214279677	;
z	=	3.280806061	;
iva	=	-2.087851335	;
aran	=	-3.565633332	;
dm	=	0.792025694	;
pi 	=	0.029558802	;
vr	=	0.036093325	;
ppet	=	0	;
qpet	=	-0.220768927	;
lamda	=	-0.051293294	;
ga	=	-0.332478425	;
po	=	0	;
v1	=	3.272312421	;
v2	=	3.323605715	;
iri	=	0.097264242	;
gb	=	-0.974937334	;
rb	=	-0.49866339	;
irb	=	0.055490357	;
tc	=	1.92449695	;
cb	=	-1.710224263	;
ira	=	0.063256546	;
walras	=	0	;
scc	=	-0.117259771	;
sck	=	0.119359758	;
pw	=	0	;
iv	=	-0.127625826	;
rv	=	-3.041199347	;
pv	=	0	;
yv	=	-0.127625826	;
yt	=	2.59992291	;
pt	=	0	;
taukv	=	0	;
kapa	=	0	;
ltv	=	-0.356674944	;
irh	=	0.113328685	;
ch	=	-0.013498939	;
tauv	=	0	;
rrb	=	-2.566081655	;
psiv	=	6.003657866	;
psim	=	-2.606793879	;
psin	=	-0.174441902	;
u	=	5.444972737	;
pomeg	=	-2.301550103	;
omeg	=	-0.049735903	;
rs	=	-2.162823151	;
sdomeg	=	-3.248639373	;
chi	=	-6.362528228	;

end;

shocks;
 %var e1 = sigma^2;
% var e2 = sigma^2;
% var e3 = sigma^2;
 %var e4 = sigma^2;
% var e5 = sigma^2;
% var e6 = sigma^2;
% var e7 = sigma^2;
% var e8 = sigma^2;
% var e9 = sigma^2;
% var e10= sigma^2;
% var e11= sigma^2;
% var e12= sigma^2;
%var e13= sigma^2;
%var e14= sigma^2;
% var e15= sigma^2;
%var e16=sigma^2;
% var e17=sigma^2;
 var e18=sigma^2;
% var e19= sigma^2;
% var e29= sigma^2;
% var e21= sigma^2;
% var e22=sigma^2;
%var e23=sigma^2;
%var e24= sigma^2;
% var e25=sigma^2;
% var e26=sigma^2;

%corr e18,e24=0.8 ;


// Se activan y se aplican las Correlaciones simultaneas de  choques aleatorios.
%var e3,e4= cor1*sigma*sigma; 
%var e3,e15= cor1*sigma*sigma; 
%var e4,e15=cor1*sigma*sigma;
%var e24,e25=corl*sigma*sigma;
%var e24,e24=cor1*sigma*sigma;
end;

%corr e1, e2=0.8;

resid(1);
steady;
check;

%stoch_simul(hp_filter = 2100, order = 1);
%stoch_simul(hp_filter =400, order = 1, irf=20) pib,pi,er,c,i,vr,irb;
%stoch_simul(hp_filter =400, order = 1, irf=20);
%stoch_simul (periods=2100) z, f, xt, vr, irb, ifdi, ppet, qpet, pib, n, pi, er, c, i, g, dm;

%stoch_simul (periods=2100, order=1)z,ppet,ifdi,xt,pib,n,pi,dm,er,c,i,g,f,vr,iva,aran,tauk,irb,qpet,pwx,fk,tauh,x,m,y,ir,iri,irpi,iripi,df,s;


stoch_simul (periods=2100, order=1);
%stoch_simul (periods=2100, order=2);
dynasave(RESULT) pib, y ;


%Otras ecuaciones que se han usado:
%  exp(kv)=(1/psiv)^(-1/gamav)*(kai*exp(c)^(-gamac)-beta*theta*kai*exp(c(+1))^(-gamac)*((exp(tauv))*exp(rv)+exp(taukv)*exp(pv)*deltav/theta))^(-1/gamav);
%  exp(cb(-1))/exp(bf(-1))=(deltatc*(exp(irb)/exp(ir)))^(-sigmatc);
%  exp(ir)=exp(irb)+spread2;
%  exp(iri(+1))+exp(rrb)/(1-exp(rrb))-exp(ira(+1))/(1-exp(rrb))=kf*alfaf*exp(df)^(alfaf-1);
%  exp(iri)=exp(ira)+spread;
% Supuesto de país pequeño (31)
%  exp(xt)=exp(x)+exp(xx);
%  exp(pwx)=exp(pw);
%  exp(cb)/exp(bf)=(deltatc*(exp(irb)/exp(ir)))^(-sigmatc);
%  exp(tc)*exp(ira(+1))=exp(bf)*exp(ir(+1))+exp(cb)*exp(irb(+1));
%  exp(dm)=(exp(c)^(gamac/gamam))*(exp(kapa(+1))*(exp(ir(+1))-1)/(exp(ir(+1))*psim))^(-1/gamam);
%  exp(kv)=(1/psiv)^(-1/gamav)*(exp(kapa)*exp(c)^(-gamac)*exp(taukv)*exp(pv)-beta*exp(kapa(+1))*exp(c(+1))^(-gamac)*((exp(tauv(+1)))*exp(rv(+1))
%             +exp(tauv(+1))*exp(pv(+1))*deltav))^(-1/gamav);
%  exp(s)*exp(er)=(((exp(er)*exp(s(-1)))*(1+exp(ri)))/(theta))-exp(tauk)*(((exp(r)*exp(k(-1)))/(theta)))-exp(tauh)*(exp(w)*exp(n)+exp(ga))
%               -(1-exp(tauv))*exp(rv)*exp(kv(-1))/theta-(exp(iva)*(exp(c)+exp(i)+exp(ifdi)*exp(er)+exp(g)))/(1+exp(iva))
%               -exp(aran)*pwm*exp(er)*exp(m)-(exp(taukv)-1)*exp(pv)*exp(iv)+exp(g)-share*exp(ppet)*exp(qpet)*exp(er)-oabc;
%  exp(omeg)=(((exp(bf)*exp(ir))+(exp(cb)*exp(irb))+rcc*exp(dm)-exp(rb)-act-(1-exp(rs))*exp(irh)*exp(ch)))/((1-exp(rs))*(exp(df)*exp(iri)));  
%  exp(iri)+exp(rrb)/(1-exp(rrb))-exp(ira)/(1-exp(rrb))=kf*alfaf*exp(df)^(alfaf-1)+exp(chi)*exp(pomeg);
%  exp(irh)+exp(rrb)/(1-exp(rrb))-exp(ira)/(1-exp(rrb))=kh*alfah*exp(ch)^(alfah-1)+exp(chi)*exp(pomeg);

