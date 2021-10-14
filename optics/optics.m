MProp[d_] := {{1, d}, {0, 1}};
MRefract[n1_, n2_, r_] := {{1, 0}, {(n1 - n2)/(n2 r), n1/n2}};
MWindow[d_, n1_, n2_] := 
  MRefract[n2, n1, \[Infinity]].MProp[d].MRefract[n1, n2, \[Infinity]];
MLens [f_] := {{1, 0}, {-1/f, 1}};
qM[q_, M_] := (q M[[1, 1]] + M[[1, 2]])/(q M[[2, 1]] + M[[2, 2]]);
wq[q_, z_, \[Lambda]_] := Sqrt[\[Lambda]/(\[Pi] Im[1/(q + z)])];
w[z_, w0_, \[Lambda]_] := 
  w0 (1 + ((z \[Lambda])/(\[Pi] w0^2))^2)^(1/2);

(*TODO: make beamPropagatePlot call beamPropagate*)
beamPropagate[qInit_, sys_, \[Lambda]_] := 
 Module[{q = qInit, qNext, lines = {}, i, element, z = 0, zNext, 
   waistfunc},
  (*return list of piecewise functions showing beam propagation*)
  For[i = 1, i <= Length[sys], i++,
   element = sys[[i]];
   qNext = qM[q, element[[1]]];
   zNext = z + element[[2]];
   waistfunc = 
    Piecewise[{{wq[qNext, x - z, \[Lambda]], z < x <= zNext}}];
   If[Length[waistfunc[[2]] > 1], waistfunc[[2]] = None];
   (*Print[waistfunc];*)
   waistfunc[[2]] = None;
   AppendTo[lines, waistfunc];
   q = qNext + element[[2]];
   z = zNext;
   ]; (*waist in [m]*)
  lines
  ]

beamPropagatePlot[qInit_, sys_, \[Lambda]_, logplot_: False] := 
 Module[{q = qInit, qNext, lines = {}, i, element, z = 0, zNext, x, 
   waistfunc},
  (*plots the piecewise beam propagation*)
  For[i = 1, i <= Length[sys], i++,
   element = sys[[i]];
   qNext = qM[q, element[[1]]];
   zNext = z + element[[2]];
   waistfunc = 
    Piecewise[{{wq[qNext, x - z, \[Lambda]], z < x <= zNext}}];
   If[Length[waistfunc[[2]] > 1], waistfunc[[2]] = None];
   AppendTo[lines, waistfunc];
   q = qNext + element[[2]];
   z = zNext;
   ]; (*waist in [m]*)
  If[logplot,
   LogPlot[Evaluate[lines], {x, 10*^-4, z}, PlotRange -> Automatic, 
    Frame -> False],
   Plot[lines, {x, 0, z}, PlotRange -> {0, All}, Frame -> False]
   ](*q after transformation M*)
  ]
USAFlpmm[group_, element_] := 
 2^(group + (element - 1)/6) // 
  N; (*US Air Force standard resolution [line pair/mm] where a line \
pair is one bright and one dark line*)