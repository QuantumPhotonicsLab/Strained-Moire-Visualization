(* ::Package:: *)

(* ::Input:: *)
(*WTe2 homobilayer Moire  pattern under strain*)


(* ::Input:: *)
(*Definiton of WTe2 Bravais lattice *)


(* ::Input:: *)
(*a:={0.3498,0.6338} (*Lattice constant of WTe2 in nm, https://materialsproject.org/materials/mp-22693/*)*)
(*beta:=90 (*Angle between unit vectors*)*)
(*delta:=0.0(*No lattice mismatch*)*)
(*unitWpos1[x_,y_]:={a[[1]]*x,a[[2]]*y}*)
(*unitWpos2[x_,y_]:={a[[1]]*(x+0.5),a[[2]](y+0.4)}*)
(*unitTepos1[x_,y_]:={a[[1]]*x,a[[2]]*(y+0.3)}*)
(*unitTepos2[x_,y_]:={a[[1]]*x,a[[2]]*(y+0.65)}*)
(*unitTepos3[x_,y_]:={a[[1]]*(x+0.5),a[[2]]*(y+0.14)}*)
(*unitTepos4[x_,y_]:={a[[1]]*(x+0.5),a[[2]]*(y+0.79)}*)
(*nextpos1[x_,y_]:={a[[1]]*(x+1),a[[2]]*(y+0.3)}*)
(*nextpos2[x_,y_]:={a[[1]]*(x+1),a[[2]]*(y+0.65)}*)
(*nextpos3[x_,y_]:={a[[1]]*(x+1),a[[2]]*(y+0.)}*)
(*nextpos4[x_,y_]:={a[[1]]*(x+1),a[[2]]*(y+1)}*)
(*nextpos5[x_,y_]:={a[[1]]*(x+0),a[[2]]*(y+1)}*)
(*nextpos6[x_,y_]:={a[[1]]*(x+0),a[[2]]*(y-0.3)}*)
(*unitCell[x_,y_]:={Gray,*)
(*Line[{unitWpos1[x,y],unitTepos1[x,y]}],Line[{unitWpos1[x,y],nextpos6[x,y]}],Line[{unitWpos1[x,y],unitTepos3[x,y]}],Line[{unitWpos2[x,y],unitTepos3[x,y]}],Line[{unitWpos2[x,y],unitTepos1[x,y]}],Line[{unitWpos2[x,y],unitTepos2[x,y]}],Line[{unitWpos2[x,y],unitTepos4[x,y]}],Line[{unitWpos2[x,y],nextpos1[x,y]}],Line[{unitWpos2[x,y],nextpos2[x,y]}],Line[{unitTepos3[x,y],nextpos3[x,y]}],Line[{unitTepos4[x,y],nextpos4[x,y]}],Line[{unitTepos4[x,y],nextpos5[x,y]}],Blue,Disk[unitWpos1[x,y],a[[1]]*0.1],Disk[unitWpos2[x,y],a[[1]]*0.1],Red,Disk[unitTepos1[x,y],a[[1]]*0.1],Disk[unitTepos2[x,y],a[[1]]*0.1],Disk[unitTepos3[x,y],a[[1]]*0.1],Disk[unitTepos4[x,y],a[[1]]*0.1]}*)
(*Graphics[Block[{unitVectA={1,0},unitVectB={Cos[90Degree],Sin[90 Degree]}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-2,2},{k,-2+Ceiling[j],2+Ceiling[j]}]],PlotRange->{{-1*a[[2]],1*a[[2]]},{-1*a[[2]],1*a[[2]]}}]*)
(**)


(* ::Input:: *)
(*Functions to calculate the moire vector of strained lattice*)


(* ::Input:: *)
(*psi[theta0_,beta_]:={theta0, theta0+beta}*)
(*c0[ec_,es_]:=(1+ec)^2-es^2*)
(*cmu[ec_,es_,mu_]:=(1+mu*ec)^2-mu^2*es^2*)
(*strainmatrix[ec_,es_,mu_,phis_]:={{(1+mu *ec)+mu*es*Sin[phis Degree],mu*es*Cos[phis Degree]},{mu*es*Cos[phis Degree],(1+mu *ec)-mu*es*Sin[phis Degree]}}*)
(*rotationmatrix[angle_]:={{Cos[angle Degree],-Sin[angle Degree]},{Sin[angle Degree],Cos[angle Degree]}}*)
(*nominator[ec_,es_,delta_,theta_,mu_]:=cmu[ec,es,mu]*(1+delta)^2+c0[ec,es]-2*(1+delta)*((1+mu*ec+ec)+mu*(ec^2-es^2))*Cos[theta Degree]*)
(*moirevector1[ec_,es_,delta_,theta_,mu_,phis_,theta0_,beta_]:=a[[1]]*(1+delta)/nominator[ec,es,delta,theta,mu]*(cmu[ec,es,mu]*(1+delta)*strainmatrix[ec,es,1,phis]-c0[ec,es]*strainmatrix[ec,es,mu,phis].rotationmatrix[theta]).rotationmatrix[psi[theta0,beta][[1]]].{1,0}*)
(*moirevector2[ec_,es_,delta_,theta_,mu_,phis_,theta0_,beta_]:=a[[2]]*(1+delta)/nominator[ec,es,delta,theta,mu]*(cmu[ec,es,mu]*(1+delta)*strainmatrix[ec,es,1,phis]-c0[ec,es]*strainmatrix[ec,es,mu,phis].rotationmatrix[theta]).rotationmatrix[psi[theta0,beta][[2]]].{1,0}*)


(* ::Input:: *)
(*Visualizing the Moire pattern in 2D with adjustable strain parameters*)


(* ::Input:: *)
(*Manipulate[With[{points=Block[{unitVectA={1,0},unitVectB={Cos[beta Degree],Sin[beta Degree]}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-2size,2size},{k,-size,size}]]},Column[{Graphics[{GeometricTransformation[points,strainmatrix[ec,es,1,phis].rotationmatrix[psi[theta0][[1]]]],GeometricTransformation[points,(1+delta)*strainmatrix[ec,es,mu,phis].rotationmatrix[psi[theta0][[1]]+theta]],FaceForm[None],EdgeForm[Red],Parallelogram[{0,0},{-moirevector1[ec,es,delta,theta,mu,phis,theta0,beta],-moirevector2[ec,es,delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{moirevector1[ec,es,delta,theta,mu,phis,theta0,beta],-moirevector2[ec,es,delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{-moirevector1[ec,es,delta,theta,mu,phis,theta0,beta],moirevector2[ec,es,delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{moirevector1[ec,es,delta,theta,mu,phis,theta0,beta],moirevector2[ec,es,delta,theta,mu,phis,theta0,beta]}]},PlotRange->{{-size*a[[2]],size*a[[2]]},{-size*a[[2]],size*a[[2]]}},ImageSize-> Large]}]],{{ec,0,"Biaxial Strain, \!\(\*SubscriptBox[\(\[Epsilon]\), \(c\)]\)"},-0.1,0.1,0.01},{{es,0,"Shear Strain, \!\(\*SubscriptBox[\(\[Epsilon]\), \(s\)]\)"},-0.1,0.1,0.01},{{theta,5,"Twist angle \[Theta]"},0,5.0},{{mu,0,"Strain transfer \[Mu]"},0,1},{{phis,0,"Shear Angle \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)"},0,90},{{theta0,0,"Lattice Orientation \!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\)"},0,360},{{size,15,"Amount of unit cells shown"},1,50}]*)


(* ::Input:: *)
(*Visualizing the Moire pattern in 2D with uniaxial strain*)


(* ::Input:: *)
(*Remove[ec,es,size,theta,theta0,mu,phis]*)
(*nu:=0.18(*Poisson ratio for WTe2 https://aip.scitation.org/doi/pdf/10.1063/1.4774090*)*)
(*phis:=90*)
(*ec[eu_]:=0.5*(1-nu)*eu*)
(*es[eu_]:=0.5*(1+nu)*eu*)
(*Manipulate[With[{points=Block[{unitVectA={1,0},unitVectB={Cos[beta Degree],Sin[beta Degree]}},Table[unitCell@@(unitVectA j+unitVectB k),{j,-2size,2size},{k,-size,size}]]},Column[{Graphics[{GeometricTransformation[points,strainmatrix[ec[eu],es[eu],1,phis].rotationmatrix[psi[theta0][[1]]]],GeometricTransformation[points,(1+delta)*strainmatrix[ec[eu],es[eu],mu,phis].rotationmatrix[psi[theta0][[1]]+theta]],FaceForm[None],EdgeForm[Red],Parallelogram[{0,0},{-moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta],-moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta],-moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{-moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta],moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta]}],Parallelogram[{0,0},{moirevector1[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta],moirevector2[ec[eu],es[eu],delta,theta,mu,phis,theta0,beta]}]},PlotRange->{{-size*a[[2]],size*a[[2]]},{-size*a[[2]],size*a[[2]]}},ImageSize-> Large]}]],{{eu,0,"Uniaxial Strain, \!\(\*SubscriptBox[\(\[Epsilon]\), \(u\)]\)"},-0.1,0.1,0.01},{{theta,5,"Twist angle \[Theta]"},0,5.0},{{mu,0,"Strain transfer \[Mu]"},0,1},{{phis,0,"Shear Angle \!\(\*SubscriptBox[\(\[Phi]\), \(s\)]\)"},0,90},{{theta0,0,"Lattice Orientation \!\(\*SubscriptBox[\(\[Theta]\), \(0\)]\)"},0,360},{{size,15,"Amount of unit cells shown"},1,50}]*)
