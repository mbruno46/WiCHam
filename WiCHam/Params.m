(* ::Package:: *)

(* ::Section:: *)
(*Constants and alpha_s*)


CF[Nc_:3]:=(Nc^2-1)/(2*Nc);
beta0[Nc_,Nf_]:=(11 Nc-2 Nf)/3;
beta1[Nc_,Nf_]:=(34 Nc^2-10 Nc Nf)/3-2 (Nc^2-1)/(2 Nc) Nf;
beta2[Nc_,Nf_]:=2857/54*Nc^3+CF[Nc]^2*Nf-205/18*CF[Nc]*Nc*Nf-1415/54*Nc^2*Nf+11/9*CF[Nc]*Nf^2+79/54*Nc*Nf^2;
beta3[Nc_,Nf_]:=Module[{zeta3},zeta3=1.202056903;(149753/6+3564*zeta3)-(1078361/162+6508/27*zeta3)*Nf+(50065/162+6472/81*zeta3)*Nf^2+1093/729*Nf^3];
alphasMZ=0.117;
MW=80.2;
MZ=80.2/Sqrt[1-0.23^2];
mtop=170;
mbottom=4.4;
mcharm=1.3;


alphas[mu_,L_,Nc_,Nf_,loop_:2]:=Module[{t,aLO,aa,b0,b1,b2,b3},t=Log[mu^2/L^2];
b0=beta0[Nc,Nf];b1=beta1[Nc,Nf];b2=beta2[Nc,Nf];b3=beta3[Nc,Nf];
aLO=(4*Pi)/(b0*t);
aa=aLO*{1,-b1*Log[t]/(b0^2*t),
(b1^2*(Log[t]^2-Log[t]-1)+b0*b2)/(b0^4*t^2),
-(b1^3*(Log[t]^3-2.5*Log[t]^2-2*Log[t]+0.5)+3*b0*b1*b2*Log[t]-0.5*b0^2*b3)/(b0^6*t^3)};
Sum[aa[[i]],{i,loop}]];


FindLambda[a_,mu_,Nc_,Nf_,loop_]:=Module[{rule,L},
rule=FindRoot[a-alphas[mu,x,Nc,Nf,loop],{x,0.25}];
L=x/.rule];
