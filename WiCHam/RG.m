(* ::Package:: *)

(* ::Section:: *)
(*Renormalization Group*)


U0[a1_,a2_,b0_,g0_]:=Module[{eig,V,Vinv,aux,x},Vinv=Eigenvectors[g0]//N;
V=Inverse[Vinv]//N;
eig=Eigenvalues[g0]//N;
eig=eig/(2*b0);
aux=Table[If[(eig[[j]]-eig[[i]]==1),heps,0],{i,Length[eig]},{j,Length[eig]}];
eig=eig+Table[1,{i,Length[eig]}].aux;
x=a1/a2;
Chop[Normal[Series[V.DiagonalMatrix[Table[x^eig[[i]],{i,Length[eig]}]].Vinv,{heps,0,1}]]]];
(*MatrixExp[Transpose[g0]*Log[a1/a2]/(2*b0)];*)
J[b0_,b1_,g0_,g1_]:=Module[{eig,V,Vinv,G,aux,H,J},Vinv=Eigenvectors[g0]//N;
V=Inverse[Vinv]//N;
G=Vinv.Transpose[g1].V;
eig=Eigenvalues[g0];
eig=eig/(2*b0);
aux=Table[If[eig[[j]]-eig[[i]]==1,1/(2*b0*heps),1/(2*b0)/(1+eig[[i]]-eig[[j]])],{i,Length[eig]},{j,Length[eig]}];
J=Transpose[g0]*b1/(2*b0^2)-V.(G*aux).Vinv;
Chop[Expand[J]]];
U[a1_,a2_,b0_,g0_,J_]:=Module[{u0,res},
u0=U0[a1,a2,b0,g0];
(*tmp=(IdentityMatrix[Dimensions[g0]]+a2/(4\[Pi])*J).U0[a1,a2,b0,g0].(IdentityMatrix[Dimensions[g0]]-a1/(4\[Pi])*J);*)
(*u0=MatrixExp[Transpose[g0]*Log[a1/a2]/(2*b0)];*)
res=(a2*J.u0-a1*u0.J)/(4*Pi);
res=Chop[Expand[res]]/.{heps->0};
(u0/.{heps->0})+NLO*res];


M1[b0_,b1_,ge0_,gse1_,J_]:=Transpose[gse1]-b1/b0*Transpose[ge0]+Transpose[ge0].J-J.Transpose[ge0];
R[a1_,a2_,b0_,g0_,ge0_,M1_,J_]:=Module[{R,K0,K11,Ks,M0,MM1,H,V,Vinv,eig,aux,x},eig=Eigenvalues[g0];
Vinv=Eigenvectors[g0]//N;
eig=eig/(2*b0);
V=Inverse[Vinv]//N;
M0=Vinv.Transpose[ge0].V//N;
x=a1/a2;
aux=Table[If[eig[[i]]-eig[[j]]==1,-x^eig[[j]]*Log[x]/a2,(x^eig[[j]]/a2-x^eig[[i]]/a1)/(eig[[i]]-eig[[j]]-1)],{i,Length[eig]},{j,Length[eig]}];
K0=M0*aux;
H=Chop[Expand[Vinv.J.V]];
MM1=Vinv.M1.V;
aux=Table[If[eig[[i]]==eig[[j]],-x^eig[[j]]*Log[x],(x^eig[[j]]-x^eig[[i]])/(eig[[i]]-eig[[j]])],{i,Length[eig]},{j,Length[eig]}];
K11=MM1*aux;
Ks=K11-a1*K0.H+a2*H.K0;
R=-(2*Pi/b0)*V.(K0+NLO*Ks/(4*Pi)).Vinv;
Chop[Expand[R]]/.{heps->0}];


FullU[a1_,a2_,ae_,b0_,g0_,ge0_,M1_,J_]:=U[a1,a2,b0,g0,J]+(ae/(4*Pi))*R[a1,a2,b0,g0,ge0,M1,J];


P={{0,0,-1/3,1,-1/3,1,0,0,0,0}};
deltarsMb=5/18*Transpose[P].{{0,0,0,-2,0,-2,0,1,0,1}};
deltarsMc=-5/9*Transpose[P].{{0,0,0,1,0,1,0,1,0,1}};
Phat={{0,0,0,0,0,0,1,0,1,0}};
deltareMb=10/81*Transpose[Phat].{{0,0,6,2,6,2,-3,-1,-3,-1}};
deltareMc=-40/81*Transpose[Phat].{{0,0,3,1,3,1,3,1,3,1}};
M[mu_,a_,ae_:0,size_:10]:=Module[{res},If[mu==mbottom,
res=IdentityMatrix[size]+NLO*deltarsMb[[1;;size,1;;size]]*(a/(4*Pi))+NLO*deltareMb[[1;;size,1;;size]]*(ae/(4*Pi));];
If[mu==mcharm,
res=IdentityMatrix[size]+NLO*deltarsMc[[1;;size,1;;size]]*(a/(4*Pi))+NLO*deltareMc[[1;;size,1;;size]]*(ae/(4*Pi));];
Return[res]];


Fs[a_,z12_]:=Module[{z1,z2},z1=z12[[1]];z2=z12[[2]];
Transpose[{{z1,z2,NLO*a/(36*Pi)*z2,-NLO*a/(12*Pi)*z2,NLO*a/(36*Pi)*z2,-NLO*a/(12*Pi)*z2}}]];
Fse[a_,z12_,aw_:0]:=Module[{z1,z2,Fe},z1=z12[[1]];z2=z12[[2]];
Fe=-4/9*(3*z1+z2);
Transpose[{{z1,z2,NLO*a/(36*Pi)*z2,-NLO*a/(12*Pi)*z2,NLO*a/(36*Pi)*z2,-NLO*a/(12*Pi)*z2,NLO*aw/(6*Pi)*Fe,0,NLO*aw/(6*Pi)*Fe,0}}]];


ReduceOrder[expr_,order_:"LO"]:=Module[{tmp,expr2},
expr2=Normal[Series[expr,{NLO,0,1}]];
If[order=="LO",expr2=expr2/.{NLO->0}];
If[order=="NLO",expr2=expr2/.{NLO->1}];
expr2=Normal[Series[expr2,{as,0,1},{ae,0,1}]];
expr2=Expand[expr2/.{ae->t}/.{as->t}];
expr2=Normal[Series[expr2,{t,0,1}]];
Return[expr2/.{t->1}]];

