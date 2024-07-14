#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*Network in the Continuum*)


(* ::Subsubsection:: *)
(*dual problem definition*)


(* ::Input:: *)
(**)
(*c=1/Sqrt[1+(Grad[Subscript[p, dual][x,y],{x,y}])^2]*)
(*(*c=List[1,1];*)*)
(*opdual=Inactive[Div][({{c[[1]],0},{0,c[[2]]}} . Inactive[Grad][Subscript[p, dual][x,y],{x,y}]),{x,y}]*)


(* ::Subsubsection:: *)
(*Geometry definition*)


(* ::Input:: *)
(*r=0.005;*)
(*in=((x+1)^2+(y+1)^2<=r^2);*)
(*ground = ((x-1)^2+(y-1)^2<=r^2);*)
(*X = ((x+1)^2+(y-1)^2<=r^2);*)
(*Y = ((x-1)^2+(y+1)^2<=r^2);*)
(*\[CapitalOmega]dual=ImplicitRegion[!*)
(*X&&!Y&&!in &&!ground,{{x,-2,2},{y,-2,2}}];*)
(*\[CapitalOmega]=ImplicitRegion[!in &&!ground,{{x,-2,2},{y,-2,2}}];*)


(* ::Input:: *)
(*\[CapitalOmega]dual*)
(*\[CapitalOmega]*)


(* ::Subsubsection:: *)
(*BCs measurement*)


(* ::Input:: *)
(*Subscript[\[CapitalGamma], N]=NeumannValue[0 ,{x==2,y==2,x==-2,y==-2}];*)
(*pDirichletConds = {DirichletCondition[p[x,y]==1,in],*)
(*Subscript[\[CapitalGamma], D,ground]=DirichletCondition[p[x,y]==0,ground]};*)
(*(*Subscript[\[CapitalGamma], D,in]=DirichletCondition[p[x,y]==1,in];*)
(*Subscript[\[CapitalGamma], D,ground]=DirichletCondition[p[x,y]==0,ground];*)*)


(* ::Subsubsection:: *)
(*Functions*)


(* ::Input:: *)
(*ErrorFunc[A_,B_,px_,py_]:={A-px,B-py}*)
(*InpNew[InpOld_,\[Alpha]_,Error_]:=InpOld-\[Alpha](Error[[1]]+Error[[2]])*)
(*xpNew[xpOld_,\[Alpha]_,Error_]:={xpOld[[1]]+\[Alpha](Error[[1]]),xpOld[[2]]+\[Alpha](Error[[2]])}*)


(* ::Subsection:: *)
(*Loop*)


(* ::Input:: *)
(*(*Iterate to solve for p_dual multiple times*)*)
(*iterations=5;*)
(*pduals=Table[0,{iterations}];*)
(*pduals[[1]]=pdual0;*)
(*(*Subscript[Inp, dual]=1.039;*)*)
(*Subscript[Inp, dual]=1.0;*)
(*(*Subscript[xp, dual]={0.2,0.23};*)*)
(*Subscript[xp, dual]={0.2,0.8};*)
(*(*Error={-0.069,0.030};*)*)
(*Error={0.0,0.0};*)
(*InpLst={Subscript[Inp, dual]};*)
(*xpLst={Subscript[xp, dual]};*)
(*ErrorLst={Error};*)
(*\[Alpha]=1.0;*)
(*A=0.25;*)
(*B=0.35;*)


(* ::Input:: *)
(*ErrorLst*)


(* ::Subsubsection:: *)
(*BCs dual*)


(* ::Input:: *)
(**)
(*Subscript[Inp, dual]=InpNew[InpLst[[-1]],\[Alpha],ErrorLst[[-1]]]*)
(*Subscript[xp, dual]=xpNew[xpLst[[-1]],\[Alpha],ErrorLst[[-1]]]*)


(* ::Input:: *)
(*pdualDirichletConds = {DirichletCondition[Subscript[p, dual][x,y]==Subscript[Inp, dual],in],*)
(*DirichletCondition[Subscript[p, dual][x,y]==0,ground],*)
(*DirichletCondition[Subscript[p, dual][x,y]==Subscript[xp, dual][[1]],X],*)
(*DirichletCondition[Subscript[p, dual][x,y]==Subscript[xp, dual][[2]],Y]};*)
(*(*Subscript[\[CapitalGamma], D,in,dual]=DirichletCondition[Subscript[p, dual][x,y]==1,in];*)
(*Subscript[\[CapitalGamma], D,ground,dual]=DirichletCondition[Subscript[p, dual][x,y]==0,ground];*)
(*Subscript[\[CapitalGamma], D,X,dual]=DirichletCondition[Subscript[p, dual][x,y]==0.3,X];*)
(*Subscript[\[CapitalGamma], D,Y,dual]=DirichletCondition[Subscript[p, dual][x,y]==0.2,Y];*)*)


(* ::Input:: *)
(*(*pdualfun=NDSolveValue[{Subscript[op, dual]==Subscript[\[CapitalGamma], N],Subscript[\[CapitalGamma], D,in,dual],Subscript[\[CapitalGamma], D,ground,dual],Subscript[\[CapitalGamma], D,X,dual],Subscript[\[CapitalGamma], D,Y,dual]},Subscript[p, dual],{x,y}\[Element]\[CapitalOmega]];*)*)
(*Clear[pfun];*)
(*Clear[op];*)
(*pdualfun=NDSolveValue[{opdual==Subscript[\[CapitalGamma], N],pdualDirichletConds},Subscript[p, dual],{x,y}\[Element]\[CapitalOmega]dual];*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Plot3D[Evaluate[pdualfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(p\), \(dual\)]\)",16,Bold],AxesLabel->{"x","y","p"}]*)
(*(*Plot3D[pfun[x,y],{x,y}\[Element]\[CapitalOmega]]*)*)


(* ::Input:: *)
(*gradPdual=Grad[pdualfun[x,y],{x,y}]*)
(*\[Sigma][x_,y_]:=1/Sqrt[1+gradPdual^2];*)
(**)
(*(*\[Sigma][x_,y_]:=1/Sqrt[1+(Grad[pdualfun[x,y],{x,y}])^2]*)*)
(*Plot3D[Evaluate[\[Sigma][x,y][[1]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(xx\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)
(*Plot3D[Evaluate[\[Sigma][x,y][[2]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(yy\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)


(* ::Input:: *)
(**)
(*pdualfun[x,y]*)
(*(*\[Sigma][x,y] /. {x->1,y->1}*)*)


(* ::Input:: *)
(*(*Define sigma as an interpolated function*)*)
(*points=Flatten[Table[{x,y},{x,-2,2,0.1},{y,-2,2,0.1}],1];*)


(* ::Input:: *)
(*1/Sqrt[1+gradPdual[[1]]/. {x->0.1,y->0.1}]*)
(*Table[{x,y,1},{x,y}\[Element]points]*)


(* ::Input:: *)
(*sigmaValues1=Table[{x,y,1/Sqrt[1+gradPdual[[1]]/. {x->x,y->y}]},{x,y}\[Element]points];*)
(*sigmaValues2=Table[{x,y,1/Sqrt[1+gradPdual[[2]]/. {x->x,y->y}]},{x,y}\[Element]points];*)
(**)


(* ::Input:: *)
(*sigmaFunc1=Interpolation[sigmaValues1];*)


(* ::Input:: *)
(**)
(**)
(*(*op=Inactive[Div][({{\[Sigma][x,y][[1]],0},{0,\[Sigma][x,y][[2]]}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)
(*op=Inactive[Div][({{sigmaFunc[x,y][[1]],0},{0,sigmaFunc[x,y][[2]]}} . Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)
(*op*)
(*opdual*)
(**)


(* ::Input:: *)
(**)
(*(*pfun=NDSolveValue[{op==Subscript[\[CapitalGamma], N],{pDirichletConds}},p,{x,y}\[Element]\[CapitalOmega]];*)*)
(*pfun=NDSolve[{op==Subscript[\[CapitalGamma], N],{pDirichletConds}},p,{x,y}\[Element]\[CapitalOmega]];*)


(* ::Input:: *)
(*(*Plot3D[Evaluate[pfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)*)


(* ::Input:: *)
(*Plot3D[Evaluate[p[x,y]/.pfun],{x,-2,2},{y,-2,2}, RegionFunction->Function[{x,y,z},RegionMember[Subscript[\[CapitalOmega], dual],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)


(* ::Input:: *)
(*Evaluate[p[0.2,0.3]/.pfun]*)


(* ::Input:: *)
(*px={pfun[-1+r,1+r],pfun[1+r,-1+r]}*)
(*Error=ErrorFunc[A,B,px[[1]],px[[2]]]*)


(* ::Input:: *)
(*AppendTo[InpLst,Subscript[Inp, dual]];*)
(*AppendTo[xpLst, Subscript[xp, dual]];*)
(*AppendTo[ErrorLst,Error];*)


(* ::Input:: *)
(*(*Reset the recursion limit if changed*)*)
(*$RecursionLimit=1024;*)
