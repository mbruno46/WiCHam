(* ::Package:: *)

(* ::Section:: *)
(*Wilson Coeffcients*)


SinThetaW2=0.23;
C1[a_]:=11/2*a/(4 \[Pi])*NLO;
C2[a_,ae_]:=1-11/6*a/(4 \[Pi])*NLO-35/18*ae/(4*Pi)*NLO;
E0[x_]:=-(2/3)*Log[x]+(x*(18-11*x-x^2))/(12*(1-x)^3)+(x^2*(15-16*x+4*x^2))/(6*(1-x)^4)*Log[x];
B0[x_]:=0.25*(x/(1-x)+x*Log[x]/(x-1)^2);
C0[x_]:=x/8*((x-6)/(x-1)+(3*x+2)*Log[x]/(x-1)^2);
C3[a_,mt_,MW_,ae_]:=-(a/(24 \[Pi]))*(E0[mt^2/MW^2]-2/3*NLO)+ae/(6*Pi*SinThetaW2)*(2*B0[mt^2/MW^2]+C0[mt^2/MW^2]);
C4[a_,mt_,MW_]:=a/(8 \[Pi])*(E0[mt^2/MW^2]-2/3*NLO);
C5[a_,mt_,MW_]:=-(a/(24 \[Pi]))*(E0[mt^2/MW^2]-2/3*NLO);
C6[a_,mt_,MW_]:=C4[a,mt,MW];
D0[x_]:=-4/9*Log[x]+(25*x^2-19*x^3)/36/(x-1)^3+x^2*Log[x]*(5*x^2-2*x-6)/18/(x-1)^4;
C7[mt_,MW_,ae_]:=ae/(6*Pi)*(4*C0[mt^2/MW^2]+D0[mt^2/MW^2]-4/9*NLO);
C8=0;
C9[mt_,MW_,ae_]:=ae/(6*Pi)*(4*C0[mt^2/MW^2]+D0[mt^2/MW^2]-4/9*NLO+(10*B0[mt^2/MW^2]-4*C0[mt^2/MW^2])/SinThetaW2);
C10=0;


(* ::Section:: *)
(*Computation of z and y coefficients*)


ComputeZ[mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init12_:{0,0}]:=Module[{aref,z12,z,Lam,JJ,MM,Umat},
Lam=FindLambda[initAlphaMZ,Mz,3,5,loop];
JJ=J[beta0[3,5],beta1[3,5],gammas0[3,5,2,2],gammas1[5,2,2]];
MM=M1[beta0[3,5],beta1[3,5],gammae0[3,5,2,2],gammase1[5,2,2],JJ];
If[init12=={0,0},z12={C1[as*initAlphaMZ],C2[as*initAlphaMZ,ae*aem]},z12=init12];
(* above bottom *)
If[mu>=mbottom,
Umat=FullU[as*initAlphaMZ,as*alphas[mu,Lam,3,5,loop],ae*aem,beta0[3,5],gammas0[3,5,2,2],gammae0[3,5,2,2],MM,JJ];
z12=Umat . z12;
Return[Fs[0,z12]]];
(* below bottom *)
aref=alphas[mbottom,Lam,3,5,loop];
Umat=FullU[as*initAlphaMZ,as*aref,ae*aem,beta0[3,5],gammas0[3,5,2,2],gammae0[3,5,2,2],MM,JJ];
z12=Umat . z12;
(* Matching 5 to 4 flavor theories*)
Lam=FindLambda[aref,mbottom,3,4,loop];
JJ=J[beta0[3,4],beta1[3,4],gammas0[3,4,2],gammas1[4,2]];
MM=M1[beta0[3,4],beta1[3,4],gammae0[3,4,2,2],gammase1[4,2,2],JJ];
(* above charm *)
If[mu>mcharm,Umat=FullU[as*aref,as*alphas[mu,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4,2],gammae0[3,4,2,2],MM,JJ];
z12=Umat . z12;Return[Fs[0,z12]]];
(* below charm *)
Umat=FullU[as*aref,as*alphas[mcharm,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4,2],gammae0[3,4,2,2],MM,JJ];
z12=Umat . z12;
aref=alphas[mcharm,Lam,3,4,loop];
z=FseMat[as*aref,ae*aem] . z12;
If[mu==mcharm,Return[z]];
(* Matching 4 to 3 flavor theories*)
Lam=FindLambda[aref,mcharm,3,3,loop];
JJ=J[beta0[3,3],beta1[3,3],gammas0[3,3,10,1],gammas1[3,10,1]];
MM=M1[beta0[3,3],beta1[3,3],gammae0[3,3,10,1],gammase1[3,10,1],JJ];
Umat=FullU[as*aref,as*alphas[mu,Lam,3,3,loop],ae*aem,beta0[3,3],gammas0[3,3,10,1],gammae0[3,3,10,1],MM,JJ];
z=Umat . z];


ComputeV[mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init_:{0,0,0,0,0,0,0,0,0,0}]:=Module[{aref,v,Lam,JJ,MM,Umat},
Lam=FindLambda[initAlphaMZ,Mz,3,5,loop];
JJ=J[beta0[3,5],beta1[3,5],gammas0[3,5],gammas1[5]];
MM=M1[beta0[3,5],beta1[3,5],gammae0[3,5],gammase1[5],JJ];
If[init=={0,0,0,0,0,0,0,0,0,0},v={C1[as*initAlphaMZ],C2[as*initAlphaMZ,ae*aem],C3[as*initAlphaMZ,mtop,Mz,ae*aem],
C4[as*initAlphaMZ,mtop,Mz],C5[as*initAlphaMZ,mtop,Mz],C6[as*initAlphaMZ,mtop,Mz],C7[mtop,Mz,ae*aem],C8,C9[mtop,Mz,ae*aem],C10},v=init];
(* above bottom *)
If[mu>mbottom,Umat=FullU[as*initAlphaMZ,as*alphas[mu,Lam,3,5,loop],ae*aem,beta0[3,5],gammas0[3,5],gammae0[3,5],MM,JJ];
Return[Umat . v]];
(* below bottom *)
aref=alphas[mbottom,Lam,3,5,loop];
Umat=FullU[as*initAlphaMZ,as*aref,ae*aem,beta0[3,5],gammas0[3,5],gammae0[3,5],MM,JJ];
v=M[mbottom,as*aref,ae*aem] . Umat . v;
If[mu==mbottom,Return[v]];
(* Matching 5 to 4 flavor theories*)
Lam=FindLambda[aref,mbottom,3,4,loop];
JJ=J[beta0[3,4],beta1[3,4],gammas0[3,4],gammas1[4]];
MM=M1[beta0[3,4],beta1[3,4],gammae0[3,4],gammase1[4],JJ];
(* above charm *)
If[mu>mcharm,Umat=FullU[as*aref,as*alphas[mu,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4],gammae0[3,4],MM,JJ];
v=Umat . v;Return[v]];
(* below charm *)
Umat=FullU[as*aref,as*alphas[mcharm,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4],gammae0[3,4],MM,JJ];
aref=alphas[mcharm,Lam,3,4,loop];
v=M[mcharm,as*aref,ae*aem] . Umat . v;
If[mu==mcharm,Return[v]];
(* Matching 4 to 3 flavor theories*)
Lam=FindLambda[aref,mcharm,3,3,loop];
JJ=J[beta0[3,3],beta1[3,3],gammas0[3,3,10,1],gammas1[3,10,1]];
MM=M1[beta0[3,3],beta1[3,3],gammae0[3,3,10,1],gammase1[3,10,1],JJ];
Umat=FullU[as*aref,as*alphas[mu,Lam,3,3,loop],ae*aem,beta0[3,3],gammas0[3,3,10,1],gammae0[3,3,10,1],MM,JJ];
v=Umat . v];


ComputeY[z_,mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init_:vuoto]:=z-ComputeV[mu,initAlphaMZ,loop,Mz,aem,init];

