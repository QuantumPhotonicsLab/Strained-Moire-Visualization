(* ::Package:: *)

(* ::Input:: *)
(*Graphene-hBN heterobilayer Moire  pattern under strain*)


(* ::Input:: *)
(*Definiton of Graphene and hBN Bravais lattice *)


(* ::Input:: *)
(*a:=0.246(*Lattice constant of Graphene in nm, https://www.nature.com/articles/nphys2954*)*)
(*delta:=0.018(*Lattice mismatch of Graphene-hBN , https://www.nature.com/articles/nphys2954*)*)
(*unitvectorW[x_,y_]:=a*{x,y}*)
(*unitvectorSe[x_,y_]:=a*{x,y+2/3 Sin[120 Degree]}*)
(*unitCell[x_,y_]:={Black,Disk[unitvectorW[x,y],0.1a],Black,Disk[unitvectorSe[x,y],0.1a],Gray,*)
(*,Line[{a*{x,y},a*{x,y+2/3 Sin[120 Degree]}}],Line[{a*{x,y},a*{x+Cos[30 Degree]/2,y-Sin[30 Degree]/2}}],Line[{a*{x,y},a*{x-Cos[30 Degree]/2,y-Sin[30 Degree]/2}}]}*)
(*Graphics[Block[{unitVectA={Cos[120 Degree],Sin[120 Degree]},unitVectB={1,0}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-1,1},{k,-1+Ceiling[j],1+Ceiling[j]}]]]*)
(**)


(* ::Input:: *)
(*Functions to calculate the moire vector of strained lattice*)


(* ::Input:: *)
(*psi[theta0_]:={theta0, theta0+60}*)
(*c0[ec_,es_]:=(1+ec)^2-es^2*)
(*cmu[ec_,es_,mu_]:=(1+mu*ec)^2-mu^2*es^2*)
(*strainmatrix[ec_,es_,mu_,phis_]:={{(1+mu *ec)+mu*es*Sin[phis Degree],mu*es*Cos[phis Degree]},{mu*es*Cos[phis Degree],(1+mu *ec)-mu*es*Sin[phis Degree]}}*)
(*rotationmatrix[angle_]:={{Cos[angle Degree],-Sin[angle Degree]},{Sin[angle Degree],Cos[angle Degree]}}*)
(*nominator[ec_,es_,delta_,theta_,mu_]:=cmu[ec,es,mu]*(1+delta)^2+c0[ec,es]-2*(1+delta)*((1+mu*ec+ec)+mu*(ec^2-es^2))*Cos[theta Degree]*)
(*moirevector1[ec_,es_,delta_,theta_,mu_,phis_,theta0_]:=a*(1+delta)/nominator[ec,es,delta,theta,mu]*(cmu[ec,es,mu]*(1+delta)*strainmatrix[ec,es,1,phis].rotationmatrix[psi[theta0][[1]]].{1,0}-c0[ec,es]*strainmatrix[ec,es,mu,phis].rotationmatrix[psi[theta0][[1]]+theta].{1,0})*)
(*moirevector2[ec_,es_,delta_,theta_,mu_,phis_,theta0_]:=a*(1+delta)/nominator[ec,es,delta,theta,mu]*(cmu[ec,es,mu]*(1+delta)*strainmatrix[ec,es,1,phis].rotationmatrix[psi[theta0][[2]]].{1,0}-c0[ec,es]*strainmatrix[ec,es,mu,phis].rotationmatrix[psi[theta0][[2]]+theta].{1,0})*)


(* ::Input:: *)
(*Visualizing the Moire pattern in 2D with adjustable strain parameters*)


(* ::Input:: *)
(*Manipulate[With[{points=Block[{unitVectA={Cos[120 Degree],Sin[120 Degree]},unitVectB={1,0}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-size,size},{k,-size+Ceiling[j],size+Ceiling[j]}]]},Column[{Graphics[{GeometricTransformation[points,strainmatrix[ec,es,1,phis].rotationmatrix[psi[theta0][[1]]]],GeometricTransformation[points,(1+delta)*strainmatrix[ec,es,mu,phis].rotationmatrix[psi[theta0][[1]]+theta]],FaceForm[None],EdgeForm[Red],Parallelogram[{0,0},{-moirevector1[ec,es,delta,theta,mu,phis,theta0],-moirevector2[ec,es,delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{moirevector1[ec,es,delta,theta,mu,phis,theta0],-moirevector2[ec,es,delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{-moirevector1[ec,es,delta,theta,mu,phis,theta0],moirevector2[ec,es,delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{moirevector1[ec,es,delta,theta,mu,phis,theta0],moirevector2[ec,es,delta,theta,mu,phis,theta0]}]},PlotRange->{{-a*1.66size,a*1.66size},{-a*size,a*size}},ImageSize-> Large]}]],{{ec,0,"Biaxial Strain, \!\(\*SubscriptBox[\(\[Epsilon]\), \(c\)]\)"},-0.5,0.5},{{es,0,"Shear Strain, \!\(\*SubscriptBox[\(\[Epsilon]\), \(s\)]\)"},-0.5,0.5},{{theta,1,"Twist angle \[Theta]"},0,5.1},{{mu,0,"Strain transfer \[Mu]"},0,1},{{phis,0,"Shear Angle \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)"},0,90},{{theta0,0,"Lattice Orientation \!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\)"},0,360},{{size,5,"Amount of unit cells shown"},1,50}]*)


(* ::Input:: *)
(*Visualizing the Moire pattern in 2D with uniaxial strain*)


(* ::Input:: *)
(*Remove[ec,es,size,theta,theta0,mu,phis]*)
(*nu:=0.19(*Poisson ratio for Graphene https://link.springer.com/article/10.1007/s12274-014-0691-9*)*)
(*phis:=90*)
(*ec[eu_]:=0.5*(1-nu)*eu*)
(*es[eu_]:=0.5*(1+nu)*eu*)
(*Manipulate[With[{points=Block[{unitVectA={Cos[120 Degree],Sin[120 Degree]},unitVectB={1,0}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-size,size},{k,-size+Ceiling[j],size+Ceiling[j]}]]},Column[{Graphics[{GeometricTransformation[points,strainmatrix[ec[eu],es[eu],1,phis].rotationmatrix[psi[theta0][[1]]]],GeometricTransformation[points,(1+delta)*strainmatrix[ec[eu],es[eu],mu,phis].rotationmatrix[psi[theta0][[1]]+theta]],FaceForm[None],EdgeForm[Red],Parallelogram[{0,0},{-moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0],-moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0],-moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{-moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0],moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0]}],Parallelogram[{0,0},{moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0],moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0]}]},PlotRange->{{-1.66*a*size,1.66*a*size},{-1*a*size,1*a*size}},ImageSize-> Large]}]],{{eu,0,"Uniaxial, \!\(\*SubscriptBox[\(\[Epsilon]\), \(u\)]\)"},-0.1,0.1,0.01},{{theta,1,"Twist angle \[Theta]"},0,5.1},{{mu,0,"Strain transfer \[Mu]"},0,1},{{phis,0,"Shear Angle \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)"},0,90},{{theta0,0,"Lattice Orientation \!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\)"},0,360},{{size,5,"Amount of unit cells shown"},1,50}]*)
