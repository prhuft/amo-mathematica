(* ::Package:: *)

INuc = 3/2;(* nuclear spin *)
Se = 1/2;(* electron spin *)
mRb = 1.4192261*10^-25; (*[kg]*)
gS = 2.00023; (*Lande g, spin*)
gL = 1; 
gI = -0.000995;
D2MatElem = 3.584*10^-29;(*<5S1/2||er||5P3/2>[C*m]*)
D1MatElem = 2.537*10^-29;(*<5S1/2||er||5P1/2>[C*m]*)
IsatD2 = 2.5033; (*effective far-detuned. mW/cm^2*)
IsatD1 = 4.4845; (*effective far-detuned. mW/cm^2*)
IsatD2SI = 25.033; (*effective far-detuned. W/m^2*)
IsatD1SI = 44.845; (*effective far-detuned. W/m^2*)
\[Lambda]D1 = 794.97885098*10^-9;(*[m]*) 
\[Lambda]D2 = 780.24120968613*10^-9;(*[m]*)
\[Nu]D1 = 377.10746354; (*[THz]*)
\[Nu]D2 = 384.230484468562;(*[THz]*)
\[Nu]HF = 6.83468261090429 ;(*[GHz]*)
\[CapitalGamma]D2 = 2 \[Pi] 6.0659*10^6 ;(*[rad Hz]*)
\[CapitalGamma]D1 = 2 \[Pi] 5.7468*10^6 ;(*[rad Hz]*)
EHF = h 10^9 \[Nu]HF;(*[J]*)
\[Nu]D2GHz=384230.484468562;
\[Nu]D1GHz=377107.46354;
Rb87fGHz[{n_,L_,J_,F_}]:=Switch[{n,L,J,F},{5,0,1/2,1},-4.27167663181519,{5,0,1/2,2},2.563005979089,{5,1,1/2,1},\[Nu]D1GHz-0.51041019,{5,1,1/2,2},\[Nu]D1GHz+0.30624611,{5,1,3/2,0},\[Nu]D2GHz-0.302073888,{5,1,3/2,1},\[Nu]D2GHz-0.229851856,{5,1,3/2,2},\[Nu]D2GHz-0.072911332,{5,1,3/2,3},\[Nu]D2GHz+.193740846,_,0];
