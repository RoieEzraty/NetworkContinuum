#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*Network in the Continuum*)


(* ::Subsubsection:: *)
(*dual problem definition*)


(* ::Input:: *)
(**)
(*c=1/Sqrt[1+(Grad[pdual[x,y],{x,y}])^2]*)
(*(*c=List[1,10];*)*)
(*opdual=Inactive[Div][({{c[[1]],0},{0,c[[2]]}} . Inactive[Grad][pdual[x,y],{x,y}]),{x,y}]*)
(*op = Inactive[Div][({{1, 0}, {0, 1}} . Inactive[Grad][p[x, y], {x, y}]), {x, y}]*)
(*(*op=Inactive[Div][({{sigmaxPoly,0},{0,sigmayPoly}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)


(* ::Subsubsection:: *)
(*Geometry definition*)


(* ::Input:: *)
(*r=0.1;*)
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
(*Neumann BCs measurement*)


(* ::Input:: *)
(*Subscript[\[CapitalGamma], N]=NeumannValue[0 ,{x==2,y==2,x==-2,y==-2}];*)


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
(*Inpdual=1.0;*)
(*(*Subscript[xp, dual]={0.2,0.23};*)*)
(*xpdual={0.2,0.8};*)
(*(*Error={-0.069,0.030};*)*)
(*Error={0.0,0.0};*)
(*InpLst={Inpdual};*)
(*xpLst={xpdual};*)
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
(*Inpdual=InpNew[InpLst[[-1]],\[Alpha],ErrorLst[[-1]]]*)
(*xpdual=xpNew[xpLst[[-1]],\[Alpha],ErrorLst[[-1]]]*)


(* ::Input:: *)
(*pdualDirichletConds = {DirichletCondition[pdual[x,y]==Inpdual,in],*)
(*DirichletCondition[pdual[x,y]==0,ground],*)
(*DirichletCondition[pdual[x,y]==xpdual[[1]],X],*)
(*DirichletCondition[pdual[x,y]==xpdual[[2]],Y]};*)
(*(*Subscript[\[CapitalGamma], D,in,dual]=DirichletCondition[Subscript[p, dual][x,y]==1,in];*)
(*Subscript[\[CapitalGamma], D,ground,dual]=DirichletCondition[Subscript[p, dual][x,y]==0,ground];*)
(*Subscript[\[CapitalGamma], D,X,dual]=DirichletCondition[Subscript[p, dual][x,y]==0.3,X];*)
(*Subscript[\[CapitalGamma], D,Y,dual]=DirichletCondition[Subscript[p, dual][x,y]==0.2,Y];*)*)


(* ::Input:: *)
(*(*pdualfun=NDSolveValue[{Subscript[op, dual]==Subscript[\[CapitalGamma], N],Subscript[\[CapitalGamma], D,in,dual],Subscript[\[CapitalGamma], D,ground,dual],Subscript[\[CapitalGamma], D,X,dual],Subscript[\[CapitalGamma], D,Y,dual]},Subscript[p, dual],{x,y}\[Element]\[CapitalOmega]];*)*)
(*(*Clear[pfun];*)
(*Clear[op];*)*)
(*pdualfun=NDSolveValue[{opdual==Subscript[\[CapitalGamma], N],pdualDirichletConds},pdual,{x,y}\[Element]\[CapitalOmega]dual]*)
(*pfun=NDSolveValue[{op==Subscript[\[CapitalGamma], N], {DirichletCondition[p[x, y] == 1, in], DirichletCondition[p[x, y] == 0, ground]}},p,{x,y}\[Element]\[CapitalOmega]];*)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Plot3D[Evaluate[pdualfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(p\), \(dual\)]\)",16,Bold],AxesLabel->{"x","y","p"}]*)
(*(*Plot3D[Evaluate[pfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)*)
(*(*Plot3D[pfun[x,y],{x,y}\[Element]\[CapitalOmega]]*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*gradPdual=Grad[pdualfun[x,y],{x,y}];*)
(*(*gradPdualx=Derivative[1,0][pdualfun]*)
(*gradPdualy=Derivative[0,1][pdualfun]*)*)
(*(*\[Sigma][x_,y_]:={1/Sqrt[1+gradPdualx^2],1/Sqrt[1+gradPdualy^2]};*)*)
(*(*\[Sigma]x = 1/Sqrt[1+Evalute[gradPdualx]^2]*)
(*\[Sigma]y=1/Sqrt[1+gradPdualy^2]*)*)
(*\[Sigma][x_,y_]:=1/Sqrt[1+(Grad[pdualfun[x,y],{x,y}])^2];*)
(*gradPdual*)
(*Plot3D[Evaluate[gradPdual[[1]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(gradPdual\), \(x\)]\)",16,Bold],AxesLabel->{"x","y","gradPdual"}]*)
(*Plot3D[Evaluate[gradPdual[[2]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(gradPdual\), \(y\)]\)",16,Bold],AxesLabel->{"x","y","gradPdual"}]*)
(**)


(* ::Input:: *)
(*(*sigmaxTable=Flatten[Table[{{x,y},1/Sqrt[1+Derivative[1,0][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)
(*sigmayTable=Flatten[Table[{{x,y},1/Sqrt[1+Derivative[0,1][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)*)
(*sigmaxTable=Flatten[Table[{x,y,1/Sqrt[1+Derivative[1,0][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)
(*sigmayTable=Flatten[Table[{x,y,1/Sqrt[1+Derivative[0,1][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)
(*(*xyTable=Flatten[Table[{x,y},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)
(*sigmaxTable=Flatten[Table[{1/Sqrt[1+Derivative[1,0][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)
(*sigmayTable=Flatten[Table[{1/Sqrt[1+Derivative[0,1][pdualfun][x,y]^2]},{x,-2,2,0.2},{y,-2,2,0.2}],1];*)*)
(*(*sigmaxTable = 1/Sqrt[1+dxpdualfunTable^2]*)
(*sigmayTable = 1/Sqrt[1+dypdualfunTable^2]*)*)


(* ::Input:: *)
(*cleansigmaxTable=Select[sigmaxTable,FreeQ[#,Indeterminate]&];*)
(*cleansigmayTable=Select[sigmayTable,FreeQ[#,Indeterminate]&];*)


(* ::Input:: *)
(*(*data={{1,33},{2,80},{5,286},{10,771},{15,1382},{20,2087},{25,2867},{30,3707},{40,5526},{50,7470},{60,9482},{70,11507},{80,13495},{90,15391},{100,17313},{110,18631},{120,19752},{125,20064}};*)
(*sigmaxPoly=LinearModelFit[data,{x^3, x^2,x,1},{x}];*)*)
(*(*data[[1]]*)*)
(*(*cleansigmaxTable[[1]]*)*)
(*(*Length[xyTable]*)*)
(*(*sigmaxPoly=NonlinearModelFit[xyTable,sigmaxTable,{a x^3, b y^3,d x^2,e y^2,f yx,g y,h x,i},{a, b, d, e, f, g, h, i},{x,y}];*)*)
(*(*sigmaxPoly=NonlinearModelFit[cleansigmaxTable,a x^3+ b y^3+d x^2+e y^2+f y x,{a, b, d, e,f},{x,y}]*)*)
(*sigmaxPoly=Normal[NonlinearModelFit[cleansigmaxTable,a/(b+(x-1)^2+(y-1)^2)+d/(e+(x+1)^2+(y-1)^2)+f/(g+(x-1)^2+(y+1)^2)+h/(i+(x+1)^2+(y+1)^2)+1,{a, b, d, e, f, g, h, i},{x,y}]]*)
(*sigmayPoly=Normal[NonlinearModelFit[cleansigmayTable,a/(b+(x-1)^2+(y-1)^2)+d/(e+(x+1)^2+(y-1)^2)+f/(g+(x-1)^2+(y+1)^2)+h/(i+(x+1)^2+(y+1)^2)+1,{a, b, d, e, f, g, h, i},{x,y}]]*)
(*(*sigmaxPoly=InterpolatingPolynomial[cleansigmaxTable,{x,y}]*)
(*sigmayPoly=InterpolatingPolynomial[cleansigmayTable,{x,y}]*)*)


(* ::Input:: *)
(*Plot3D[sigmaxPoly,{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(xx\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)
(*Plot3D[sigmayPoly,{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(yy\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)
(*Plot3D[Evaluate[\[Sigma][x,y][[1]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(xx\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)
(*Plot3D[Evaluate[\[Sigma][x,y][[2]]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega]dual,{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["\!\(\*SubscriptBox[\(\[Sigma]\), \(yy\)]\)",16,Bold],AxesLabel->{"x","y","\[Sigma]"}]*)


(* ::Subsubsection:: *)
(*Dirichlet BC for normal problem*)


(* ::Input:: *)
(*pDirichletConds = {DirichletCondition[p[x, y] == 1, in], DirichletCondition[p[x, y] == 0, ground]};*)
(*(*op=Inactive[Div][({{Evaluate[\[Sigma][x,y][[1]]],0},{0,Evaluate[\[Sigma][x,y][[2]]]}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)
(*(*op=Inactive[Div][({{sigmaFunc[x,y][[1]],0},{0,sigmaFunc[x,y][[2]]}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)
(**)


(* ::Subsubsection:: *)
(*Equation normal problem*)


(* ::Input:: *)
(*op = Inactive[Div][({{1, 0}, {0, 1}} . Inactive[Grad][p[x, y], {x, y}]), {x, y}];*)
(*(*op=Inactive[Div][({{sigmaxPoly,0},{0,sigmayPoly}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)
(*op*)


(* ::Subsubsection:: *)
(*Solve normal problem*)


(* ::Input:: *)
(*p*)
(*Subscript[\[CapitalGamma], N]*)
(*in*)
(*ground*)
(*\[CapitalOmega]*)
(*x*)
(*y*)
(*pfun=NDSolveValue[{op==Subscript[\[CapitalGamma], N], {DirichletCondition[p[x, y] == 1, in], DirichletCondition[p[x, y] == 0, ground]}},p,{x,y}\[Element]\[CapitalOmega]];*)


(* ::Subsubsection:: *)
(*Plot normal problem*)


(* ::Input:: *)
(*Plot3D[Evaluate[p[x,y]/.pfun],{x,-2,2},{y,-2,2}, RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)


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
