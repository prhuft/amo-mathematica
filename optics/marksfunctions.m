(* ::Package:: *)

BeginPackage["MarksFunctions`","PhysicalConstants`"];

(* last edits 2010.june .10 *)

(*
MarksFunctions::usage = "define whittaker and some clebsch gordan functions."
*)


(* quantum defect energy levels *)
(* nClear[nzz,lzz,jzz]
EnljSI[nzz_,lzz_,jzz_]:=-ERySI/(nzz-qdnlfs[nzz,lzz,jzz])^2;
*)


(*
D12 is restricted Wigner rotation matrix for j=1/2 and polar rotations about z axis
d12 is reduced Wigner rotation matrix for j=1/2
*)

Clear[D12,d12,a,b,m,mp]

d12[m_,mp_,b_]:=(
 KroneckerDelta[m,1/2]KroneckerDelta[mp,1/2]Cos[b/2]
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-1/2](-1)Sin[b/2] 
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,1/2]Sin[b/2] 
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-1/2]Cos[b/2] 
 );

D12[m_,mp_,b_]:=Exp[I (Pi/2)(m-mp)]d12[m,mp,b];


(*
D32 is restricted Wigner rotation matrix for j=3/2 and polar rotations about z axis

d32 is reduced Wigner rotation matrix for j=3/2
*)

Clear[D32,d32,a,b,c,m,mp]

d32[m_,mp_,b_]:=(
 KroneckerDelta[m,3/2]KroneckerDelta[mp,3/2]Cos[b/2]^3
+KroneckerDelta[m,3/2]KroneckerDelta[mp,1/2](-1)Sqrt[3]Sin[b/2]Cos[b/2]^2 
+KroneckerDelta[m,3/2]KroneckerDelta[mp,-1/2]Sqrt[3]Sin[b/2]^2 Cos[b/2] 
+KroneckerDelta[m,3/2]KroneckerDelta[mp,-3/2](-1)Sin[b/2]^3 

+ KroneckerDelta[m,1/2]KroneckerDelta[mp,3/2]Sqrt[3]Sin[b/2]Cos[b/2]^2 
+KroneckerDelta[m,1/2]KroneckerDelta[mp,1/2](3 Cos[b/2]^2-2)Cos[b/2] 
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-1/2](3 Sin[b/2]^2-2)Sin[b/2] 
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-3/2]Sqrt[3]Sin[b/2]^2 Cos[b/2] 

 + KroneckerDelta[m,-1/2]KroneckerDelta[mp,3/2]Sqrt[3] Sin[b/2]^2 Cos[b/2]
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,1/2](-1)(3 Sin[b/2]^2-2)Sin[b/2] 
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-1/2](3 Cos[b/2]^2-2)Cos[b/2]  
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-3/2](-1)Sqrt[3]Sin[b/2]Cos[b/2]^2 

+ KroneckerDelta[m,-3/2]KroneckerDelta[mp,3/2]Sin[b/2]^3
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,1/2]Sqrt[3]Sin[b/2]^2 Cos[b/2] 
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,-1/2]Sqrt[3]Sin[b/2] Cos[b/2]^2 
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,-3/2]Cos[b/2]^3 );

D32[m_,mp_,b_]:=Exp[I (Pi/2)(m-mp)]d32[m,mp,b];


Clear[D52,d52,a,b,c,m,mp]

d52[m_,mp_,b_]:=( KroneckerDelta[m,5/2]KroneckerDelta[mp,5/2]Cos[b/2]^5
+KroneckerDelta[m,5/2]KroneckerDelta[mp,3/2](-1)Sqrt[5]Sin[b/2]Cos[b/2]^4
+KroneckerDelta[m,5/2]KroneckerDelta[mp,1/2]Sqrt[10]Sin[b/2]^2Cos[b/2]^3 
+KroneckerDelta[m,5/2]KroneckerDelta[mp,-1/2](-1)Sqrt[10]Sin[b/2]^3 Cos[b/2]^2 
+KroneckerDelta[m,5/2]KroneckerDelta[mp,-3/2]Sqrt[5]Sin[b/2]^4Cos[b/2] 
+KroneckerDelta[m,5/2]KroneckerDelta[mp,-5/2](-1)Sin[b/2]^5 

+KroneckerDelta[m,3/2]KroneckerDelta[mp,5/2] Sqrt[5]Sin[b/2]Cos[b/2]^4
+KroneckerDelta[m,3/2]KroneckerDelta[mp,3/2]Cos[b/2]^3(1-5 Sin[b/2]^2)
+KroneckerDelta[m,3/2]KroneckerDelta[mp,1/2](-1)Sqrt[2]Sin[b/2]Cos[b/2]^2(2-5 Sin[b/2]^2) 
+KroneckerDelta[m,3/2]KroneckerDelta[mp,-1/2](-1)Sqrt[2]Sin[b/2]^2 Cos[b/2](2-5 Cos[b/2]^2)
+KroneckerDelta[m,3/2]KroneckerDelta[mp,-3/2]Sin[b/2]^3 (1-5 Cos[b/2]^2)
+KroneckerDelta[m,3/2]KroneckerDelta[mp,-5/2] Sqrt[5]Sin[b/2]^4Cos[b/2] 

+KroneckerDelta[m,1/2]KroneckerDelta[mp,5/2]Sqrt[10]Sin[b/2]^2Cos[b/2]^3 
+KroneckerDelta[m,1/2]KroneckerDelta[mp,3/2]Sqrt[2]Sin[b/2]Cos[b/2]^2(2-5 Sin[b/2]^2)
+KroneckerDelta[m,1/2]KroneckerDelta[mp,1/2]Cos[b/2](3-12 Cos[b/2]^2+10 Cos[b/2]^4)           
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-1/2](-1)Sin[b/2](3-12 Sin[b/2]^2+10 Sin[b/2]^4)
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-3/2](-1)Sqrt[2]Sin[b/2]^2 Cos[b/2](2-5 Cos[b/2]^2)
+KroneckerDelta[m,1/2]KroneckerDelta[mp,-5/2](-1)Sqrt[10]Sin[b/2]^3 Cos[b/2]^2 

+KroneckerDelta[m,-1/2]KroneckerDelta[mp,5/2]Sqrt[10]Sin[b/2]^3 Cos[b/2]^2 
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,3/2](-1)Sqrt[2]Sin[b/2]^2 Cos[b/2](2-5 Cos[b/2]^2)
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,1/2]Sin[b/2](3-12 Sin[b/2]^2+10 Sin[b/2]^4)
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-1/2]Cos[b/2](3-12 Cos[b/2]^2+10 Cos[b/2]^4) 
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-3/2] (-1)Sqrt[2]Sin[b/2]Cos[b/2]^2(2-5 Sin[b/2]^2)
+KroneckerDelta[m,-1/2]KroneckerDelta[mp,-5/2]Sqrt[10]Sin[b/2]^2 Cos[b/2]^3

+KroneckerDelta[m,-3/2]KroneckerDelta[mp,5/2]Sqrt[5]Sin[b/2]^4Cos[b/2] 
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,3/2](-1)Sin[b/2]^3 (1-5 Cos[b/2]^2)
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,1/2](-1)Sqrt[2]Sin[b/2]^2 Cos[b/2](2-5 Cos[b/2]^2)
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,-1/2] Sqrt[2]Sin[b/2]Cos[b/2]^2(2-5 Sin[b/2]^2) 
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,-3/2] Cos[b/2]^3(1-5 Sin[b/2]^2)
+KroneckerDelta[m,-3/2]KroneckerDelta[mp,-5/2](-1)Sqrt[5]Sin[b/2]Cos[b/2]^4

+KroneckerDelta[m,-5/2]KroneckerDelta[mp,5/2]Sin[b/2]^5 
+KroneckerDelta[m,-5/2]KroneckerDelta[mp,3/2]Sqrt[5]Sin[b/2]^4Cos[b/2]
+KroneckerDelta[m,-5/2]KroneckerDelta[mp,1/2] Sqrt[10]Sin[b/2]^3 Cos[b/2]^2 
+KroneckerDelta[m,-5/2]KroneckerDelta[mp,-1/2]Sqrt[10]Sin[b/2]^2 Cos[b/2]^3 
+KroneckerDelta[m,-5/2]KroneckerDelta[mp,-3/2]Sqrt[5]Sin[b/2]Cos[b/2]^4
+KroneckerDelta[m,-5/2]KroneckerDelta[mp,-5/2]Cos[b/2]^5 );

D52[m_,mp_,b_]:=Exp[I (Pi/2)(m-mp)]d52[m,mp,b];



(*
Define normalized radial wave functions
Normalization from Seaton, Mon. Not. R. Astron. Soc. 188, 504 (1969)
\[Nu] is effective quantum number (\[Nu]=1/E^2 with E in Rydbergs)
l is orbital angular momentum quantum number
r is distance in atomic units
*)


Clear[kappa,mu,z,r]
Clear[nzzz,\[Alpha]zzz,xzzz]
MyLaguerreL[nzzz_,\[Alpha]zzz_,xzzz_]:=Module[{L0,L1,LL,j},L1=0.;LL=1.;Do[L0=L1;L1=LL;LL=Expand[((2j-1.+\[Alpha]zzz-xzzz)L1-(j-1.+\[Alpha]zzz)L0)/j];,{j,1,nzzz}];LL]
(* MyHypU valid for nzzz<-1 *)
MyHypU[nzzz_,\[Alpha]zzz_,xzzz_]:=Module[{index,L0,L1,LL,j},index=-Floor[-nzzz];
L1=HypergeometricU[nzzz-index,\[Alpha]zzz,xzzz];
LL=HypergeometricU[nzzz-index-1,\[Alpha]zzz,xzzz];
Do[L0=L1;L1=LL;LL=(2(nzzz-index-j)-\[Alpha]zzz+xzzz+2)L1-((nzzz-index-j)+1)((nzzz-index-j)-\[Alpha]zzz+2)L0;,{j,2,-index}];LL];

wh[kappazz_, muzz_, z_]:=Exp[-z/2] z^(muzz + 1/2) HypergeometricU[1/2 + muzz - kappazz, 1 + 2 muzz, z];
whN[kappazz_, muzz_, z_]:=Exp[-z/2] z^(muzz + 1/2) MyHypU[1/2 + muzz - kappazz, 1 + 2 muzz, z];

(* following form valid for integer arguments only *) 
whint[\[Kappa]zz_, \[Mu]zz_, z_]:= Exp[-z/2]z^(\[Mu]zz + 1/2) (-1)^(\[Kappa]zz - \[Mu]zz - 1/2)Factorial[\[Kappa]zz - \[Mu]zz - 1/2]LaguerreL[\[Kappa]zz - \[Mu]zz - 1/2, 2  \[Mu]zz, z];
whintN[\[Kappa]zz_, \[Mu]zz_, z_]:= Exp[-z/2]z^(\[Mu]zz + 1/2) (-1)^(\[Kappa]zz - \[Mu]zz - 1/2)Factorial[\[Kappa]zz - \[Mu]zz - 1/2]MyLaguerreL[\[Kappa]zz - \[Mu]zz - 1/2, 2  \[Mu]zz, z];
    
pnl[\[Nu]zz_, lzz_,rzz_]:= (1/Sqrt[Gamma[\[Nu]zz + lzz + 1] Gamma[\[Nu]zz - lzz] \[Nu]zz^2]) wh[\[Nu]zz, lzz + 1/2, 2 rzz/\[Nu]zz];
pnlN[\[Nu]zz_, lzz_,rzz_]:= (1/Sqrt[Gamma[\[Nu]zz + lzz + 1] Gamma[\[Nu]zz - lzz] \[Nu]zz^2]) whN[\[Nu]zz, lzz + 1/2, 2 rzz/\[Nu]zz];
pnlint[\[Nu]zz_, lzz_,rzz_]:= (1/Sqrt[Gamma[\[Nu]zz + lzz + 1] Gamma[\[Nu]zz - lzz] \[Nu]zz^2]) whint[\[Nu]zz, lzz + 1/2, 2 rzz/\[Nu]zz];
pnlintN[\[Nu]zz_, lzz_,rzz_]:= (1/Sqrt[Gamma[\[Nu]zz + lzz + 1] Gamma[\[Nu]zz - lzz] \[Nu]zz^2]) whintN[\[Nu]zz, lzz + 1/2, 2 rzz/\[Nu]zz];
      
   
(* **************************************************** *)
(* hyperfine transition factor squared *)
Clear[jezz, lezz, fezz, mezz, jgzz, lgzz, fgzz, mgzz, Inuczz, szz, qzz];

rsqhyp[fezz_, mezz_, fgzz_, mgzz_, qzz_, jezz_, lezz_, jgzz_, lgzz_, Inuczz_,szz_] := ((2fgzz + 1)(2fezz + 1)(2jgzz + 1)(2jezz + 1)Max[lgzz,lezz]ThreeJSymbol[{fezz, -mezz}, {1, qzz}, {fgzz, mgzz}]^2 SixJSymbol[{jezz, fezz, Inuczz}, {fgzz, jgzz, 1}]^2 SixJSymbol[{lezz, jezz, szz}, {jgzz,lgzz, 1}]^2);

chf1to2[q_,f1_,mf1_,j1_,l1_,f2_,mf2_,j2_,l2_,in_]:=((-1)^(f2-mf2)ThreeJSymbol[{f2,-mf2},{1,q},{f1,mf1}](-1)^(1+in+j2+f1) Sqrt[(2f1+1)(2f2+1)]SixJSymbol[{j2,f2,in},{f1,j1,1}](-1)^(3/2+l2+j1)Sqrt[(2j1+1)(2j2+1)]SixJSymbol[{l2,j2,1/2},{j1,l1,1}](-1)^(l2+(1+l1+l2)/2)Sqrt[Max[l1,l2]]);

mehfs12[q_,f1_,mf1_,j1_,l1_,f2_,mf2_,j2_,l2_,s_,in_]:=ClebschGordan[{f1,mf1},{1,q},{f2,mf2}](-1)^(1+in+j2+f1) Sqrt[(2f1+1)]SixJSymbol[{j1,in,f1},{f2,1,j2}](-1)^(1+s+l2+j1) Sqrt[(2j1+1)(2j2+1)]SixJSymbol[{l1,s,j1},{j2,1,l2}](-1)^(l2+(1+l1+l2)/2)Sqrt[Max[l1,l2]];  
mefs12[q_,j1_,mj1_,l1_,j2_,mj2_,l2_,s_]:=ClebschGordan[{j1,mj1},{1,q},{j2,mj2}](-1)^(1+s+l2+j1) Sqrt[(2j1+1)]SixJSymbol[{l1,s,j1},{j2,1,l2}](-1)^(l2+(1+l1+l2)/2)Sqrt[Max[l1,l2]];  

 
 
cfs1to2[q_,j1_,mj1_,l1_,j2_,mj2_,l2_]:=((-1)^(j2-mj2)ThreeJSymbol[{j2,-mj2},{1,q},{j1,mj1}](-1)^(3/2+l2+j1) Sqrt[(2j1+1)(2j2+1)]SixJSymbol[{l2,j2,1/2},{j1,l1,1}](-1)^(l2+(1+l1+l2)/2)Sqrt[Max[l1,l2]]);

(* same as previous function but reordered input arguments *)
cfs1to2b[q_,j1_,l1_,mj1_,j2_,l2_,mj2_]:=((-1)^(j2-mj2)ThreeJSymbol[{j2,-mj2},{1,q},{j1,mj1}](-1)^(3/2+l2+j1) Sqrt[(2j1+1)(2j2+1)]SixJSymbol[{l2,j2,1/2},{j1,l1,1}](-1)^(l2+(1+l1+l2)/2)Sqrt[Max[l1,l2]]);


(* coefficient of fine structure to fine structure transition with radial matrix element in J basis *)
cfs1to2J[q_,j1_,l1_,mj1_,j2_,l2_,mj2_]:=(ClebschGordan[{j1,mj1},{1,q},{j2,mj2}]/ Sqrt[2j2+1]);
 
(* --------------- Plasma dispersion function *)
Clear[z,Zpd]        
Zpd[z_]:=I Sqrt[\[Pi]] Exp[-z^2] Erfc[-I z]

(* ABCD matrices *)

Clear[z,nz,q,f,r,\[Lambda],t,n1,n2,L,\[Alpha],R];
mlinear[z_]:={{1,z},{0,1}};
mlens[f_]:={{1,0},{-1/f,1}};
mmirrortran[r_,t_,n1_,n2_]:={{1,t n1 /n2},{-(1/r)(1-n2/n1),1-(t n1/(n2 r))(1-n2/n1)}};
(* transmission through plano concave mirror of thickness t and index n2 in medium with index n1 propagating from plano to concave side  *)



mflat[n1_,n2_]:={{1,0},{0,n1/n2}};   (* transmission from index n1 to index n2 *)
(* transmission through surface of radius R and tilt angle \[Alpha] *)
mflatt[n1_,n2_,R_,\[Alpha]_]:={{Sqrt[(n2/n1)^2- Sin[\[Alpha]]^2]/((n2/n1)Cos[\[Alpha]]),0},
{(Cos[\[Alpha]]-Sqrt[n2^2/n1^2 - Sin[\[Alpha]]^2])/(R Cos[\[Alpha]] Sqrt[n2^2/n1^2 - Sin[\[Alpha]]^2]),Cos[\[Alpha]]/Sqrt[(n2/n1)^2-Sin[\[Alpha]]^2]}}
mflats[n1_,n2_,R_,\[Alpha]_]:={{1,0},
{(Cos[\[Alpha]]-Sqrt[n2^2/n1^2 - Sin[\[Alpha]]^2])/(R  n2/n1),n1/n2}}

(* reflection from curved mirror *)
mmirror[R_]:={{1,0},{-2/R,1}};
mmirrort[R_,\[Alpha]_]:={{1,0},{-2 /(Cos[\[Alpha]]R),1}}
mmirrors[R_,\[Alpha]_]:={{1,0},{-2  Cos[\[Alpha] ]/R,1}}

(* Brewster plate, thickness L *)
mBt[n1_,n2_,L_]:={{1,L n1^3/n2^3},{0,1}};
mBs[n1_,n2_,L_]:={{1,L n1/n2},{0,1}};

MyConjugate[z_]:=z/.{I->-I};
MyRe[z_]:= (1/2)(z+MyConjugate[z]);
MyIm[z_]:=(-I/2)(z-MyConjugate[z]);
wofq[q_,\[Lambda]_,nz_]:=Simplify[ Sqrt[\[Lambda]/(Pi nz Im[ComplexExpand[1/q]])],{\[Lambda]>0,nz>0,\[Lambda] \[Element] Reals, nz \[Element] Reals}];
Rofq[q_]:=1/Re[ComplexExpand[1/q]];


EndPackage[];

