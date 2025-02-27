#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*Network in the Continuum*)


(* ::Subsection:: *)
(*Preliminaries*)


(* ::Subsubsection:: *)
(*Equation definitions*)


(* ::Input:: *)
(*(*initiate conductivity field*)*)
(*R = 1;*)
(*c=1;*)
(*(*Poisson equation*)*)
(*op=Inactive[Div][({{c,0},{0,c}} . Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)
(*(*Define Rdot using active gradient and dot product*)*)
(*Rdot[\[Gamma]_,Q_,\[Beta]_,R_]:=\[Gamma] Q[x,y]-\[Beta] (Grad[R[x,y],{x,y}]\[CenterDot]Grad[R[x,y],{x,y}])*)


(* ::Subsubsection:: *)
(*Geometry definition*)


(* ::Input:: *)
(*(*Radius around inputs and outputs*)*)
(*r=0.1;*)
(*(*input domain - exclude from geometry*)*)
(*in=((x+1)^2+(y+1)^2<=r^2);*)
(*(*output domain - exclude from geometry*)*)
(*out = ((x-1)^2+(y-1)^2<=r^2);*)
(*(*geometry at which Poisson eqn. is solved*)*)
(*\[CapitalOmega]=ImplicitRegion[!in &&!out,{{x,-2,2},{y,-2,2}}]*)


(* ::Subsubsection:: *)
(*Neumann BCs at geometry boundaries*)


(* ::Input:: *)
(*Subscript[\[CapitalGamma], N]=NeumannValue[0 ,{x==2,y==2,x==-2,y==-2}];*)


(* ::Subsubsection:: *)
(*Dirichlet BCs at input and output*)


(* ::Input:: *)
(*DirichletConds = {DirichletCondition[p[x,y]==Inval,in],*)
(*DirichletCondition[p[x,y]==Outval,out]};*)


(* ::Subsubsection:: *)
(*For polynomial fits*)


(* ::Input:: *)
(*xMin=-2;*)
(*xMax=2;*)
(*yMin=-2;*)
(*yMax=2;*)
(*n=100;*)
(*xVals=Subdivide[xMin,xMax,n];*)
(*yVals=Subdivide[yMin,yMax,n];*)
(*(*Exclude points inside the small circles from interpolation data*)*)
(*dataFilter[{x_,y_,v_}]:=!((x+1)^2+(y+1)^2<=r^2||(x-1)^2+(y-1)^2<=r^2);*)


(* ::Subsubsection:: *)
(*Scalars in the problem*)


(* ::Input:: *)
(*iterations=5;  (*# iterations in dynamics of R*)*)
(*\[Beta]=0.01; (*damping coefficient of gradients in resistance *)*)
(*\[Gamma]=4.0;  (*change in R a.f.o current coefficient*)*)
(*Inval=1.0;  (*Value of pressure at input*)*)
(*Outval=0.0;  (*Value of pressure at output*)*)
(*dt=0.1; (*Time step for evolution of conductivity c*)*)
(*basis=Flatten@Table[x^i y^j,{i,0,10},{j,0,10-i}];(*Basis for a polynomial fit*)*)
(**)


(* ::Subsection:: *)
(*Dynamics*)


(* ::Input:: *)
(*(*Iterate to solve for p_dual multiple times*)*)


(* ::Input:: *)
(*pfun=NDSolveValue[{op==Subscript[\[CapitalGamma], N], DirichletConds},p,{x,y}\[Element]\[CapitalOmega]];*)
(*Plot3D[Evaluate[pfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)


(* ::Subsubsection:: *)
(*Interpolating*)


(* ::Input:: *)
(*(*Collect discrete data of pressure}*)*)
(*(*datap=Flatten[Table[{x,y,pfun[x,y]},{x,xVals},{y,yVals}],1];*)*)
(*(*datac = Flatten[Table[{x,y,c},{x,xVals},{y,yVals}],1];*)*)
(*(*datap=Select[Flatten[Table[{x,y,pfun[x,y]},{x,xVals},{y,yVals}],1],dataFilter];*)
(*datac=Select[Flatten[Table[{x,y,c},{x,xVals},{y,yVals}],1],dataFilter];*)*)
(*(*Fitp=Fit[datap,basis,{x,y}];*)*)
(*(*Fitp=NonlinearModelFit[datap,a1 Exp[-((x-1)^2+(y-1)^2)/b1]+a2 Exp[-((x+1)^2+(y+1)^2)/b2]+c1*x+c2*y+d1,{a1,b1,a2,b2, c1, c2, d1},{x,y}]*)*)
(*(*\[Epsilon]=10^-3;*)
(*Fitp=NonlinearModelFit[datap,a1 (1/((x-1)^2+(y-1)^2+\[Epsilon]))^(1/2)+a2(1/((x+1)^2+(y+1)^2+\[Epsilon]))^(1/2)+c1*x+c2*y+d1,{a1,a2, c1, c2, d1},{x,y}];*)*)
(*(*Fitp=NonlinearModelFit[datap,a1 Exp[-(Abs[x-1]+Abs[y-1])/b1]+a2 Exp[-(Abs[x+1]+Abs[y+1])/b2]+c1,{a1,b1,a2,b2, c1},{x,y}]*)*)
(*(*Fitc=Fit[datac,basis,{x,y}]*)*)
(*(*Fitc=NonlinearModelFit[datac,a1 Exp[-((x-1)^2+(y-1)^2)/b1]+a2 Exp[-((x+1)^2+(y+1)^2)/b2]+c1*x+c2*y+d1,{a1,b1,a2,b2, c1, c2, d1},{x,y}]*)*)
(*(*Fitc=NonlinearModelFit[datac,a1 (1/((x-1)^2+(y-1)^2+\[Epsilon]))^(1/2)+a2(1/((x+1)^2+(y+1)^2+\[Epsilon]))^(1/2)+c1*x+c2*y+d1,{a1,a2, c1, c2, d1},{x,y}];*)*)
(*(*Fitc=NonlinearModelFit[datac,a1 Exp[-(Abs[x-1]+Abs[y-1])/b1]+a2 Exp[-(Abs[x+1]+Abs[y+1])/b2]+c1,{a1,b1,a2,b2, c1},{x,y}]*)*)
(**)
(*(*Plot3D[Evaluate[Normal[Fitp]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p fit",16,Bold],AxesLabel->{"x","y","p polynomial fit"}]*)*)
(**)
(*(*Compute partial derivatives:*)*)
(*(*dpdx=D[Normal[Fitp],x];*)
(*dpdy=D[Normal[Fitp],y];*)
(*dcdx = D[Normal[Fitc],x]*)
(*dcdy=D[Normal[Fitc],y]*)*)
(**)
(*dpdx=D[pfun[x,y],x];*)
(*dpdy=D[pfun[x,y],y];*)
(*(*dcdx = D[c[x,y],x];*)
(*dcdy=D[c[x,y],y];*)*)
(*(*Step A:Store the symbolic expression in a variable,say absQExpr*)*)
(*(*absQExpr=Normal[Fitc]*Expand[(dpdx^2+dpdy^2)^0.5];*)*)
(*(*absQExpr=Expand[((Normal[Fitc]*dpdx)^2+(Normal[Fitc]*dpdy)^2)^0.5];*)*)
(*absQExpr=((c*dpdx)^2+(c*dpdy)^2)^0.5*)
(*(*absQExpr=((Fitc*dpdx)^2+(Fitc*dpdy)^2)^0.5*)*)
(**)
(*Plot3D[Evaluate[absQExpr],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["||Q||",16,Bold],AxesLabel->{"x","y","||Q||"}]*)
(**)
(*(*absGradc=Expand[(dcdx^2+dcdy^2)^0.5];*)*)
(*absGradc=(dcdx^2+dcdy^2)^0.5;*)
(**)
(*Plot3D[Evaluate[absGradc],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["||Grad(c)||",16,Bold],AxesLabel->{"x","y","||Grad(c)||"}]*)
(**)
(*(*cOld=Normal[Fitc]*)*)
(*(*cOld=c;*)*)
(*cOld = c*)
(*(*c=cOld +( \[Gamma] absQExpr -\[Beta] absGradc)*dt;*)*)
(*(*c=1+ \[Gamma]* pfun*dt*)*)
(*c= cOld +\[Gamma] absQExpr*)
(**)
(*datac=Select[Flatten[Table[{x,y,c},{x,xVals},{y,yVals}],1],dataFilter];*)
(*(*Fitc=NonlinearModelFit[datac,a1 (1/((x-1)^2+(y-1)^2+\[Epsilon]))^(1/2)+a2(1/((x+1)^2+(y+1)^2+\[Epsilon]))^(1/2)+c1*x+c2*y+d1,{a1,a2, c1, c2, d1},{x,y}]*)*)
(*Fitc=Fit[datac,basis,{x,y}]*)
(**)
(*Plot3D[Evaluate[c],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["c",16,Bold],AxesLabel->{"x","y","c"}]*)
(**)
(*Plot3D[Evaluate[Fitc],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["Fit c",16,Bold],AxesLabel->{"x","y","Fit c"}]*)
(**)
(*(*op=Inactive[Div][({{c,0},{0,c}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}];*)*)
(*op=Div[({{Normal[Fitc],0},{0,Normal[Fitc]}} . Grad[p[x,y],{x,y}]),{x,y}];*)
(*pfun=NDSolveValue[{op==Subscript[\[CapitalGamma], N], DirichletConds},p,{x,y}\[Element]\[CapitalOmega]];*)
(**)
(*Plot3D[Evaluate[pfun[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)
(**)


(* ::Input:: *)
(*(*(*Collect data in the form {x,y,QMagSquared[x,y]}*)*)
(*data=Flatten[Table[{x,y,pfun1[x,y]},{x,xVals},{y,yVals}],1];*)
(*basis=Flatten@Table[x^i y^j,{i,0,4},{j,0,4-i}];*)
(*polyFitp=Fit[data,basis,{x,y}];*)
(*Plot3D[Evaluate[polyFitp],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p polynomial fit",16,Bold],AxesLabel->{"x","y","p polynomial fit"}]*)
(**)
(*(*polyFitpFunc [x_,y_]:=polyFitp *)*)
(*(*1) Compute partial derivatives:*)*)
(*dpdx=D[polyFitp,x];*)
(*dpdy=D[polyFitp,y];*)
(**)
(*(*Step A:Store the symbolic expression in a variable,say absQExpr*)*)
(*absQExpr=Expand[dpdx^2+dpdy^2]*)
(*Plot3D[Evaluate[absQExpr],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["||Q||",16,Bold],AxesLabel->{"x","y","||Q||"}]*)
(*(*op=Inactive[Div][({{R[x,y],0},{0,R[x,y]}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}]*)*)
(*op2=Inactive[Div][({{absQExpr,0},{0,absQExpr}}.Inactive[Grad][p[x,y],{x,y}]),{x,y}]*)
(*pfun2=NDSolveValue[{op2==Subscript[\[CapitalGamma], N], DirichletConds},p,{x,y}\[Element]\[CapitalOmega]];*)
(*Plot3D[Evaluate[pfun1[x,y]],{x,-2,2},{y,-2,2},RegionFunction->Function[{x,y,z},RegionMember[\[CapitalOmega],{x,y}]],PlotRange->All,Mesh->None,ColorFunction->"TemperatureMap", PlotLabel->Style["p",16,Bold],AxesLabel->{"x","y","p"}]*)*)
