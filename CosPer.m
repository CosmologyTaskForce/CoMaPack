(* ::Package:: *)

(********1*********2*********3*********4*********5*********6*********7*****
 * CosPer.m = General Relativity, Einstein, Cosmology, Cosmological Perturbation 
 * Written/Modified by MA[family name] Lei: leima137@gmail.com
 *                                   ~~Department of Phsics, Fudan University~~
 *                                   Department of Phsics and Astronomy, University of New Mexico
 *                                   http://openmetric.org
 *                                   http://neutrino.xyz
 * 
 * Background part is improved based on:
 *
 * GREAT.m = General Relativity, Einstein & All That 4 Mathematica
 * written by Tristan Hubsch: thubsch@howard.edu
 *                                    Howard University Physics
 *                                    string.howard.edu/~tristan/
 * Based on the package
 * "EinsteinTensor.m"
 * by Pekka Janhunen: pjanhune@finsun.csc.fi
 *    Finnish Meteorological Institute Geophysics Dept.
 **************************************************************************)

BeginPackage["COSPER`"]

(** List of functions and their usage. If you need more explanation, please check the examples. **)

IMetric::usage = "IMetric[g], with g an n.n-matrix (with two lower indices, i.e., \!\(\*SubscriptBox[\(g\), \(\[Alpha]\[Beta]\)]\)),
  returns the inverse metric (with two upper indices, i.e., \!\(\*SuperscriptBox[\(g\), \(\[Alpha]\[Beta]\)]\))."

Christoffel::usage = "Christoffel[g,x], with g a n.n-matrix and x
  n-vector of coordinates, gives the Christoffel symbol of the 2nd
  kind (1st upper, two lower indices, i.e., \!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Beta]\[Gamma]\), \(\[Alpha]\)]\)=\!\(\*FractionBox[\(1\), \(2\)]\)\!\(\*SuperscriptBox[\(g\), \(\[Sigma]\[Rho]\)]\)(\!\(\*SubscriptBox[\(g\), \(\[Rho]\[Mu], \[Nu]\)]\)+\!\(\*SubscriptBox[\(g\), \(\[Nu]\[Rho], \[Mu]\)]\)-\!\(\*SubscriptBox[\(g\), \(\[Mu]\[Nu], \[Rho]\)]\)))." 

Riemann::usage = "Riemann[g,x], with g a n.n-matrix and x n-vector of
  coordinates, gives the Riemann tensor (1st upper, three lower
  indices, i.e., \!\(\*SubsuperscriptBox[\(R\), \(\[Mu]\[Nu]\[Sigma]\), \(\[Rho]\)]\)=\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Mu]\[Sigma], \[Nu]\), \(\[Rho]\)]\)-\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Nu]\[Sigma], \[Mu]\), \(\[Rho]\)]\)+\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Mu]\), \(\[Lambda]\)]\)\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Nu]\[Lambda]\), \(\[Rho]\)]\)-\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Sigma]\[Nu]\), \(\[Lambda]\)]\)\!\(\*SubsuperscriptBox[\(\[CapitalGamma]\), \(\[Mu]\[Lambda]\), \(\[Rho]\)]\))."

Ricci::usage = "Ricci[g,x], with g a n.n-matrix and x n-vector of
  coordinates, gives the Ricci tensor (two lower symmetric indices)."

SCurvature::usage = "SCurvature[g,x], with g a n.n-matrix and x
  n-vector of coordinates, gives the Scalar Curvature."

EinsteinTensor::usage = "EinsteinTensor[g,x] with g a nxn-matrix
  (the metric with lower indices) and x n-vector (the coordinates)
  gives the Einstein tensor (a nxn-matrix) with lower indices."

SqRicci::usage = "SqRicci[g,x], with g a n.n-matrix and x
  n-vector of coordinates, gives the norm-square of the Ricci tensor."

SqRiemann::usage = "SqRiemann[g,x], with g a n.n-matrix and x
  n-vector of coordinates, gives the norm-square of the Riemann tensor."

HubbleL::usage = "HubbleL[\[CapitalOmega]r0,\[CapitalOmega]m0,\[CapitalOmega]x0,H0,a], with \[CapitalOmega]r0, \[CapitalOmega]m0, \[CapitalOmega]x0
  the fraction of composites in the universe and H0 the Hubble constant at a=1, gives the Hubble function."

deltaChristoffel::usage = "deltaChristoffel[g,h,x], with g a background n.n-matrix metric,
  h perturbed part of metric (n.n-matrix) and x n-vector of coordinates, gives the perturbed Christoffel symbol of the 2nd
  kind (1st upper, two lower indices)."

delChristoffel::usage = "delChristoffel[g,h,x], with g a background n.n-matrix metric,
  h perturbed part of metric (n.n-matrix) and x n-vector of coordinates, gives the perturbed Christoffel symbol of the 2nd
  kind (1st upper, two lower indices)."

deltaRiemann::usage = "deltaRiemann[g,h,x], with g a background n.n-matrix metric,
  h a perturbed n.n-matrix metric and x n-vector of coordinates, gives the perturbed Riemann tensor (1st upper, three lower
  indices)."

deltaRicci::usage = "deltaRicci[g,h,x], with g a background n.n-matrix metric,
  h the perturbed part of n.n-matrix metric and x n-vector of coordinates, gives the perturbed Ricci tensor (two lower indices/Liang or Wald)."

delRicci::usage = "delRicci[g,h,x], with g a background n.n-matrix metric,
  h the perturbed part of n.n-matrix metric and x n-vector of coordinates, gives the perturbed Ricci tensor (two lower indices/Liang or Wald)."

deltaSCurvature::usage = "deltaSCurvature[g,h,x], with g a n.n-matrix metric, h a n.n-matrix perturbed part of the metric and x
  n-vector of coordinates, gives the perturbed Scalar Curvature."

deltaEinsteinTensor::usage = "deltaEinsteinTensor[g,h,x], with g a nxn-matrix
  (the metric with lower indices), h perturbed part of metric g and x n-vector (the coordinates)
  gives the perturbed Einstein tensor (a nxn-matrix) with lower indices."


helpCOSPER::usage = "COSPER functions are: IMetric, Christoffel,
  Riemann, Ricci, SCurvature, EinsteinTensor, SqRicci, SqRiemann, HubbleL, deltaChristoffel, deltaRiemann, deltaRicci, deltaSCurvature, deltaEinsteinTensor."





Begin["`Private`"]

IMetric[metric_] := Simplify[Inverse[metric]];

Christoffel[metric_,x_]:=
  Block[{Dim, iMet, PreChristoffel, Christoffel, i, j, k},
         Dim = Length[x];
         iMet = IMetric[metric];
           (* Metric with upper indices *)
         PreChristoffel =
           Table[ D[metric[[k,i]],x[[j]]]
                + D[metric[[j,k]],x[[i]]]
                - D[metric[[i,j]],x[[k]]],
          	   {k,Dim}, {i,Dim}, {j,Dim} ];
          	 (* The \{k,ij\} Christoffel symbols *)
         PreChristoffel = Simplify[PreChristoffel];
            (* The \Gamma^k_{ij} Christoffel symbols *)
         Christoffel = (1/2) iMet . PreChristoffel;
         (* Return the Christoffel symbol: *)
         Simplify[Christoffel] ]

Riemann[metric_,x_]:=
  Block[ {Dim, iMet, ChrisSymbol, Riemann, PreRiemann,
          a, b, c, i, j, k},
          Dim = Length[x];
          iMet = IMetric[metric];
          ChrisSymbol = Christoffel[metric,x];
          PreRiemann = 
             Table[ D[ChrisSymbol[[a,i,c]],x[[b]]]
                  + Sum[ChrisSymbol[[k,i,c]]
                        * ChrisSymbol[[a,k,b]],
                        {k,Dim} ],
                    {a,Dim}, {i,Dim}, {b,Dim}, {c,Dim} ];
           	(* Riemann tensor (1st upper, three lower indices)
           	   is antisymmetrized PreRiemann: *)
          Riemann = Table[ PreRiemann[[a,i,b,c]]
                         - PreRiemann[[a,i,c,b]],
                           {a,Dim}, {i,Dim}, {b,Dim}, {c,Dim} ];
          (* Return the Riemann tensor: *)
          Simplify[Riemann] ]

Ricci[metric_,x_]:=
  Block[ {Dim, Riem, Ricci, a, i, j},
          Dim = Length[x];
          Riem = Riemann[metric,x];
          Ricci = Table[ Sum[Riem[[a,i,a,j]], {a,Dim}],
                         {i,Dim}, {j,Dim} ];
          (* Return the Ricci tensor (two lower indices): *)
          Simplify[Ricci] ]

SCurvature[metric_,x_]:=
  Block[ {Dim,iMet,CurvatureScalar,i,j},
          Dim=Length[x];
          iMet=IMetric[metric];
          CurvatureScalar=Tr[iMet.Ricci[metric,x]];
    (*Return Scalar Curvature:*)Simplify[CurvatureScalar]]

EinsteinTensor[metric_,x_]:=
  Simplify[Ricci[metric,x] - (1/2) SCurvature[metric,x] metric];

SqRicci[metric_,x_]:=
  Block[ {Dim, iMet, Ric, RRicci, i, j, k, l},
          Dim = Length[x];
          iMet = IMetric[metric];
          Ric = iMet.Ricci[metric, x];
          RRicci = Tr[Ric.Ric];
          (*Return norm - square of Ricci tensor :*)Simplify[RRicci]]

SqRiemann[metric_, x_] := Block[{Dim, iMet, Riem,
   RRiem, i, j, k, l, m, n, p, q}, Dim = Length[x];
    iMet = IMetric[metric];
    Riem = Riemann[metric, x].iMet;
    RRiem = Sum[Riem[[i, j, l, k]]*
    Riem[[l, k, i, j]], {i, Dim}, {j, Dim}, {k, Dim}, {l, Dim}];
    (*Return norm - square of Riemann tensor :*)Simplify[RRiem]]

HubbleL[\[CapitalOmega]r0_,\[CapitalOmega]m0_,\[CapitalOmega]x0_,H0_,a_] := H0 Sqrt[\[CapitalOmega]r0 a^-4+\[CapitalOmega]m0 a^-3+\[CapitalOmega]x0];

(*
Inversed metric: Overscript[g, ~]^\[Sigma]\[Rho]
*)
deltaChristoffel[metric_,deltametric_,x_]:=
Block[{Dim, iMet, ddeltametric, PredelChristoffel, temp, delChristoffel, i, j, k, l},
         Dim = Length[x];
         iMet = IMetric[metric];
           (* Metric with upper indices *)
		ddeltametric = Table[D[deltametric[[i,j]],x[[k]]],{i,Dim},{j,Dim},{k,Dim}];
		temp = Table[ddeltametric[[i,j,k]] - Sum[deltametric[[i,l]]*Christoffel[metric,x][[l,j,k]],{l,Dim}]
				 - Sum[deltametric[[l,j]]*Christoffel[metric,x][[l,i,k]],{l,Dim}],{i,Dim},{j,Dim},{k,Dim}];
				(* Subscript[\[Delta]g, ij,k]-Subscript[\[Delta]g, il] Subsuperscript[\[CapitalGamma], jk, l]-Subscript[\[Delta]g, lj] Subsuperscript[\[CapitalGamma], ik, l] *)
        PredelChristoffel = 
					Table[temp[[l,k,j]] + temp[[l,j,k]] - temp[[j,k,l]],
          	  {l,Dim},{k,Dim},{j,Dim} ];
          	 (* Subscript[\[Delta]g, ij,k]-Subscript[\[Delta]g, il] Subsuperscript[\[CapitalGamma], jk, l]-Subscript[\[Delta]g, lj] Subsuperscript[\[CapitalGamma], ik, l] *)
		 PredelChristoffel = Simplify[PredelChristoffel];
         delChristoffel = (1/2)iMet.PredelChristoffel;
         (* Dot: means g^il (temp[[l,k,j]]+temp[[l,j,k]]-temp[[j,k,l]]) *)
         Simplify[delChristoffel] ]

(*
(*!!Used for testing. This is not the right perturbed Christoffel.*)
delChristoffel[metric_, deltametric_, x_]:=
		Block[{Dim, iMet, idelMet, PredelChristoffel, delChristoffel, i, j, k},
         Dim = Length[x];
		iMet = IMetric[metric];
		idelMet = IMetric[deltametric];
		PredelChristoffel = Table[D[deltametric[[k,i]],x[[j]]] + D[deltametric[[j,k]],x[[i]]]-D[deltametric[[i,j]],x[[k]]],
							{k,Dim},{i,Dim},{j,Dim}]
		delChristoffel = 1/2 idelMet.Christoffel[metric,x]+1/2 iMet.PredelChristoffel;
         Simplify[delChristoffel] ]
*)

(*
Subscript[R, \[Mu]\[Nu]]=Subscript[\[Delta]\[CapitalGamma]^\[Nu], \[Mu]\[Sigma],\[Nu]]-Subscript[\[Delta]\[CapitalGamma]^\[Nu], \[Nu]\[Sigma],\[Mu]]+Subscript[\[Delta]\[CapitalGamma]^\[Lambda], \[Mu]\[Sigma]] Subscript[\[CapitalGamma]^\[Nu], \[Lambda]\[Nu]]+Subscript[\[CapitalGamma]^\[Lambda], \[Mu]\[Sigma]] Subscript[\[Delta]\[CapitalGamma]^\[Nu], \[Lambda]\[Nu]]-Subscript[\[Delta]\[CapitalGamma]^\[Lambda], \[Nu]\[Sigma]] Subscript[\[CapitalGamma]^\[Nu], \[Lambda]\[Mu]]-Subscript[\[CapitalGamma]^\[Lambda], \[Nu]\[Sigma]] Subscript[\[Delta]\[CapitalGamma]^\[Nu], \[Lambda]\[Mu]]
*)
deltaRicci[metric_,deltametric_,x_]:=
  Block[ {Dim, p1, p2, Gamma, deltaGamma, dGamma, deltaRicci1, deltaRicci2, deltaRicci, i, j, k, l, m, n},
          Dim = Length[x];
		deltaGamma = deltaChristoffel[metric,deltametric,x];
		Gamma = Christoffel[metric,x];
		(*dGamma is the gradient of Christoffel.*)
		deltaRicci = Table[Sum[Sum[ - deltaGamma[[l,i,k]]*Gamma[[k,j,l]] - deltaGamma[[k,j,l]]*Gamma[[l,i,k]],{l,Dim}],{k,Dim}]
					 + Sum[ - D[deltaGamma[[k,i,k]],x[[j]]] + D[deltaGamma[[k,i,j]],x[[k]]],{k,Dim}]
					 + Sum[Sum[deltaGamma[[l,i,j]]*Gamma[[k,k,l]] + deltaGamma[[k,k,l]]*Gamma[[l,i,j]],{k,Dim}],{l,Dim}],{i,Dim},{j,Dim}];
          (* Return the perturbed Ricci tensor (two lower indices): *)
(*BACKUP: (*This is another faster calculation but may give the wrong answer.*)
        deltaRicci1 = Tr[Transpose[dGamma,{1,3,2,4}],Plus,2] - Tr[Transpose[dGamma,{1,4,2,3}],Plus,2];
		deltaRicci2 = Table[Sum[deltaGamma[[l,i,k]]*Gamma[[k,j,l]],{k,Dim},{l,Dim}]
					 + Sum[deltaGamma[[k,j,l]]*Gamma[[l,i,k]],{k,Dim},{l,Dim}]
					 - Sum[deltaGamma[[l,i,j]]*Tr[Gamma,Plus,2][[l]],{l,Dim}]
					 - Sum[Tr[deltaGamma,Plus,2][[l]]*Gamma[[l,i,j]],{l,Dim}],{i,Dim},{j,Dim}];
*)          Simplify[deltaRicci] ]

(*
(*Not varified. You can adapt this for use.*)
deltaRiemann[metric_,deltametric_,x_]:=
  Block[ {Dim, iMet, ChrisSymbol, deltaGamma, deltaRiemann, PredeltaRiemann,
          a, b, c, i, j, k},
          Dim = Length[x];
          iMet = IMetric[metric];
          ChrisSymbol = Christoffel[metric,x];
		deltaGamma = deltaChristoffel[metric,deltametric,x];
          PredeltaRiemann = 
             Table[ D[deltaGamma[[a,i,c]],x[[b]]]
                  + Sum[ChrisSymbol[[k,i,c]]
                        * deltaGamma[[a,k,b]],
                        {k,Dim} ]
				 + Sum[deltaGamma[[k,i,c]]
                        * ChrisSymbol[[a,k,b]],
                        {k,Dim} ],
                    {a,Dim}, {i,Dim}, {b,Dim}, {c,Dim} ];
           	(* Perturbed Riemann tensor (1st upper, three lower indices)
           	   is antisymmetrized PreRiemann: *)
          deltaRiemann = Table[ PredeltaRiemann[[a,i,b,c]]
                         - PredeltaRiemann[[a,i,c,b]],
                           {a,Dim}, {i,Dim}, {b,Dim}, {c,Dim} ];
          (* Return the Riemann tensor: *)
          Simplify[deltaRiemann] ]
*)

(*(*Used for testing. Not the right perturbed Ricci.*)
delRicci[metric_,deltametric_,x_]:=
  Block[ {Dim, delRiem, delRicci, a, i, j},
          Dim = Length[x];
          delRiem = deltaRiemann[metric,deltametric,x];
          delRicci = Table[ Sum[delRiem[[a,i,a,j]], {a,Dim}],
                         {i,Dim}, {j,Dim} ];
          (* Return the Ricci tensor (two lower indices): *)
          Simplify[delRicci] ]
*)

deltaSCurvature[metric_,deltametric_,x_]:=
	Block[ {Dim,iMet,deltaCurvatureScalar,i,j},
          Dim=Length[x];
          iMet=IMetric[metric];
		idelMet=IMetric[deltametric];
          deltaCurvatureScalar=Sum[idelMet[[i,j]]*Ricci[metric,x][[i,j]]+iMet[[i,j]]*deltaRicci[metric,deltametric,x][[i,j]],{i,Dim},{j,Dim}];
    (*Return perturbed Scalar Curvature:*)Simplify[deltaCurvatureScalar]
		]

deltaEinsteinTensor[metric_,deltametric_,x_]:=
  Simplify[deltaRicci[metric,deltametric,x] - (1/2) SCurvature[metric,x] . deltametric
			 - (1/2) deltaSCurvature[metric,deltametric,x]. metric];


T0WMAP=2.726;
hWMAP=0.719;
N\[Nu]WMAP=3.046; (* WMAP recommends 3.04... who knows why? This value is taken from Lesgourgues&Pastor *)
\[CapitalOmega]m0h2WMAP=(0.1326);
\[CapitalOmega]b0h2WMAP=0.02273;

As2WMAP=2.41; 
nsWMAP=0.963;
\[Tau]reiWMAP=0.087;
ReionizationFractionWMAP=1.;

CPL:={T0WMAP,hWMAP,N\[Nu]WMAP,\[CapitalOmega]m0h2WMAP,\[CapitalOmega]b0h2WMAP,As2WMAP,nsWMAP,\[Tau]reiWMAP,ReionizationFractionWMAP,{},0,1};

Digit[num_]:=NumberForm[num,{5}];
Cosmology[cpl_,order_]:=Grid[
Join[{{"Variable","Value","Units","Comment"},{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(b0\)]\)",Digit[\[CapitalOmega]b0[cpl]],"","Abundance of baryons"},{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(c0\)]\)",Digit[\[CapitalOmega]c0[cpl]],"","Abundance of CDM"},
{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r0\)]\)",Digit[\[CapitalOmega]r0[cpl]],"","Abundance of radiation (massless \[Nu]'s and photons)"}},If[MassiveNeutrinosBool[cpl],{{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Nu]m0\)]\)",Digit[\[CapitalOmega]\[Nu]m0[cpl]],"","Abundance of massive neutrinos"},
{"\!\(\*SubscriptBox[\(m\), \(\[Nu]\)]\)'s",MassesNeutrinoseV[cpl],"eV","Masses of massive neutrinos"}},{}],
{{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[CapitalLambda]0\)]\)",Digit[\[CapitalOmega]\[CapitalLambda]0[cpl]],"","Abundance of \[CapitalLambda]"},{"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(K\)]\)",0,"","Abundance of curvature"},{"\!\(\*SubscriptBox[\(T\), \(0\)]\)",Digit[T0[cpl]],"K","Temperature of CMB"},{"\!\(\*SubscriptBox[\(N\), \(\[Nu]\)]\)",Digit[Nneutrinos[cpl]],"","Number of massless neutrinos"},{"h",Digit[hred[cpl]],"","Reduced Hubble constant"},{"\!\(\*SubscriptBox[\(\[Tau]\), \(rei\)]\)",Digit[\[Tau]rei[cpl]],"","Optical depth of reionization"},{"\!\(\*SubscriptBox[\(n\), \(s\)]\)",Digit[SpectralIndex[cpl]],"","Scalar perturbations spectral index"},{"\!\(\*SubscriptBox[\(k\), \(eq\)]\)",Digit[1/RescaleMp[cpl]],"\!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)","k at equivalence time"},{"\!\(\*SubscriptBox[\(z\), \(rei\)]\)",Digit[zOFy[cpl][yBrei[cpl]]],"","Redshift at reionization"},{"\!\(\*SubscriptBox[\(z\), \(eq\)]\)",Digit[zOFy[cpl][1]],"","Redshift at equivalence"},{"\!\(\*SubscriptBox[\(z\), \(LSS\)]\)",Digit[zLSS[cpl]],"","Redshift at \[Tau]-\!\(\*SubscriptBox[\(\[Tau]\), \(rei\)]\)=ln(2)"},{"\!\(\*SubscriptBox[\(z\), \(dec\)]\)",Digit[zdec[cpl]],"","Redshift at max of visibility function (\[Tau]'\!\(\*SuperscriptBox[\(e\), \(-\[Tau]\)]\)) "},{"\!\(\*SubscriptBox[\(z\), \(*\)]\)",Digit[zstar[cpl]],"","Redshift at \[Tau]-\!\(\*SubscriptBox[\(\[Tau]\), \(rei\)]\)=1"},{"\!\(\*SubscriptBox[\(d\), \(A\)]\)(\!\(\*SubscriptBox[\(z\), \(*\)]\))",Digit[Dist[cpl][etaOFz[cpl][zstar[cpl]]]*RescaleMp[cpl]],"Mpc","Angular distance at \!\(\*SubscriptBox[\(z\), \(*\)]\)"},
{"\!\(\*SubscriptBox[\(d\), \(A\)]\)(\!\(\*SubscriptBox[\(z\), \(eq\)]\))",Digit[Dist[cpl][etaOFy[cpl][1]]*RescaleMp[cpl]],"Mpc","Angular distance at equivalence"},
{"\!\(\*SubscriptBox[\(D\), \(H\)]\)",1/H[cpl][y0[cpl]]*RescaleMp[cpl],"Mpc","Hubble distance today"},{"\!\(\*SubscriptBox[\(t\), \(0\)]\)",10^(-9)AgeUniverse[cpl][y0[cpl]],"Gyears","Age of the Universe"},{"\!\(\*SubscriptBox[\(t\), \(*\)]\)",AgeUniverse[cpl][yOFz[cpl][zstar[cpl]]],"years","Age of universe at \!\(\*SubscriptBox[\(z\), \(*\)]\)"},{"\!\(\*SubscriptBox[\(r\), \(hor\)]\)(\!\(\*SubscriptBox[\(z\), \(dec\)]\))",Digit[etaOFz[cpl][zdec[cpl]]*RescaleMp[cpl]],"Mpc","Radius of horizon at \!\(\*SubscriptBox[\(z\), \(dec\)]\)"},{"\!\(\*SubscriptBox[\(\[Eta]\), \(0\)]\)",Digit[\[Eta]0[cpl]*RescaleMp[cpl]],"Mpc","Conformal time today"}},Switch[order,PerturbationParameters,{{"\!\(\*SuperscriptBox[SubscriptBox[\(A\), \(s\)], \(2\)]\)",Digit[SpectrumAmplitude[cpl]],"","Primordial scalar perturbations amplitude at k=0.002 Mpc"},
{"\!\(\*SubscriptBox[\(n\), \(S\)]\)",Digit[SpectralIndex[cpl]],"","Scalar spectral index"},{"r",Digit[rTtoS[cpl]],"","Tensor to Scalar ratio at k=0.002 Mpc"},
{"\!\(\*SubscriptBox[\(n\), \(T\)]\)",Digit[TensorSpectralIndex[cpl]],"","Tensor spectral index"},{"\!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\)",Digit[\[Sigma]8[cpl]],"","Relies on an extrapolation of the matter power spectrum if \!\(\*SubscriptBox[\(k\), \(max\)]\)<200\!\(\*SubscriptBox[\(k\), \(eq\)]\)"}},BackgroundParameters,{}]],Frame->All]


helpCOSPER:= Print["COSPER functions are: IMetric, Christoffel,
  Riemann, Ricci, SCurvature, EinsteinTensor, SqRicci, SqRiemann, HubbleL, deltaChristoffel, deltaRiemann, deltaRicci, deltaSCurvature, deltaEinsteinTensor."]


End[]

EndPackage[]

helpCOSPER

Print["Enter 'helpCOSPER' for this list of functions"]














