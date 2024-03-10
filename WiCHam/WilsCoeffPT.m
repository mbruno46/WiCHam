(* ::Package:: *)

(* ::Section:: *)
(*Computation of z and y coefficients*)


ComputeZ[mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init12_:vuoto]:=Module[{aref,z12,z,Lam,JJ,MM,Umat},
Lam=FindLambda[initAlphaMZ,Mz,3,5,loop];
JJ=J[beta0[3,5],beta1[3,5],gammas0[3,5,2,2],gammas1[5,2,2]];
MM=M1[beta0[3,5],beta1[3,5],gammae0[3,5,2,2],gammase1[5,2,2],JJ];
If[init12==vuoto,init12={C1[as*initAlphaMZ],C2[as*initAlphaMZ,ae*aem]}];
(* above bottom *)
If[mu>=mbottom,Umat=FullU[as*initAlphaMZ,as*alphas[mu,Lam,3,5,loop],ae*aem,beta0[3,5],gammas0[3,5,2,2],gammae0[3,5,2,2],MM,JJ];
Print[Umat];
z12=Umat.init12;Return[Fs[0,z12]]];
(* below bottom *)
aref=alphas[mbottom,Lam,3,5,loop];
Umat=FullU[as*initAlphaMZ,as*aref,ae*aem,beta0[3,5],gammas0[3,5,2,2],gammae0[3,5,2,2],MM,JJ];
z12=Umat.init12;
(* Matching 5 to 4 flavor theories*)
Lam=FindLambda[aref,mbottom,3,4,loop];
JJ=J[beta0[3,4],beta1[3,4],gammas0[3,4,2],gammas1[4,2]];
MM=M1[beta0[3,4],beta1[3,4],gammae0[3,4,2,2],gammase1[4,2,2],JJ];
(* above charm *)
If[mu>mcharm,Umat=FullU[as*aref,as*alphas[mu,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4,2],gammae0[3,4,2,2],MM,JJ];
z12=Umat.z12;Return[Fs[0,z12]]];
(* below charm *)
Umat=FullU[as*aref,as*alphas[mcharm,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4,2],gammae0[3,4,2,2],MM,JJ];
z12=Umat.z12;
aref=alphas[mcharm,Lam,3,4,loop];
z=Fse[as*aref,z12,ae*aem];
If[mu==mcharm,Return[z]];
(* Matching 4 to 3 flavor theories*)
Lam=FindLambda[aref,mcharm,3,3,loop];
JJ=J[beta0[3,3],beta1[3,3],gammas0[3,3,10,1],gammas1[3,10,1]];
MM=M1[beta0[3,3],beta1[3,3],gammae0[3,3,10,1],gammase1[3,10,1],JJ];
Umat=FullU[as*aref,as*alphas[mu,Lam,3,3,loop],ae*aem,beta0[3,3],gammas0[3,3,10,1],gammae0[3,3,10,1],MM,JJ];
z=Umat.z];


ComputeV[mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init_:vuoto]:=Module[{aref,v,Lam,JJ,MM,Umat},
Lam=FindLambda[initAlphaMZ,Mz,3,5,loop];
JJ=J[beta0[3,5],beta1[3,5],gammas0[3,5],gammas1[5]];
MM=M1[beta0[3,5],beta1[3,5],gammae0[3,5],gammase1[5],JJ];
If[init==vuoto,init={C1[as*initAlphaMZ],C2[as*initAlphaMZ,ae*aem],C3[as*initAlphaMZ,mtop,Mz,ae*aem],
C4[as*initAlphaMZ,mtop,Mz],C5[as*initAlphaMZ,mtop,Mz],C6[as*initAlphaMZ,mtop,Mz],C7[mtop,Mz,ae*aem],C8,C9[mtop,Mz,ae*aem],C10}];
(* above bottom *)
If[mu>mbottom,Umat=FullU[as*initAlphaMZ,as*alphas[mu,Lam,3,5,loop],ae*aem,beta0[3,5],gammas0[3,5],gammae0[3,5],MM,JJ];
Return[Umat.init]];
(* below bottom *)
aref=alphas[mbottom,Lam,3,5,loop];
Umat=FullU[as*initAlphaMZ,as*aref,ae*aem,beta0[3,5],gammas0[3,5],gammae0[3,5],MM,JJ];
v=M[mbottom,as*aref,ae*aem].Umat.init;
If[mu==mbottom,Return[v]];
(* Matching 5 to 4 flavor theories*)
Lam=FindLambda[aref,mbottom,3,4,loop];
JJ=J[beta0[3,4],beta1[3,4],gammas0[3,4],gammas1[4]];
MM=M1[beta0[3,4],beta1[3,4],gammae0[3,4],gammase1[4],JJ];
(* above charm *)
If[mu>mcharm,Umat=FullU[as*aref,as*alphas[mu,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4],gammae0[3,4],MM,JJ];
v=Umat.v;Return[v]];
(* below charm *)
Umat=FullU[as*aref,as*alphas[mcharm,Lam,3,4,loop],ae*aem,beta0[3,4],gammas0[3,4],gammae0[3,4],MM,JJ];
aref=alphas[mcharm,Lam,3,4,loop];
v=M[mcharm,as*aref,ae*aem].Umat.v;
If[mu==mcharm,Return[v]];
(* Matching 4 to 3 flavor theories*)
Lam=FindLambda[aref,mcharm,3,3,loop];
JJ=J[beta0[3,3],beta1[3,3],gammas0[3,3,10,1],gammas1[3,10,1]];
MM=M1[beta0[3,3],beta1[3,3],gammae0[3,3,10,1],gammase1[3,10,1],JJ];
Umat=FullU[as*aref,as*alphas[mu,Lam,3,3,loop],ae*aem,beta0[3,3],gammas0[3,3,10,1],gammae0[3,3,10,1],MM,JJ];
v=Umat.v];


ComputeY[z_,mu_,initAlphaMZ_,loop_,Mz_:MZ,aem_:(1/129),init_:vuoto]:=z-ComputeV[mu,initAlphaMZ,loop,Mz,aem,init];

