(* ::Package:: *)

BeginPackage[ "oscillLib`"]

(* first order solution *)

matrix[1]={{0,1,0,0},{(I \[Omega])/\[Epsilon]^2+zp/\[Epsilon]^2,0,-(zp/\[Epsilon]^2),0},{0,0,0,1},{-((zp \[Kappa])/\[Epsilon]^2),0,(I \[Omega])/(\[Epsilon]^2 \[Chi])+(zp \[Kappa])/\[Epsilon]^2,0}};
matrixEigensystem[1]=Eigensystem[matrix[1]];
expMat[1]=MatrixExp[matrix[1]t];
\[Lambda]s[1]=FullSimplify[matrixEigensystem[1][[1]][[1]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
\[Lambda]s[2]=FullSimplify[matrixEigensystem[1][[1]][[2]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
\[Lambda]s[3]=FullSimplify[matrixEigensystem[1][[1]][[3]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
\[Lambda]s[4]=FullSimplify[matrixEigensystem[1][[1]][[4]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
vs[1]=FullSimplify[matrixEigensystem[1][[2]][[1]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
vs[2]=FullSimplify[matrixEigensystem[1][[2]][[2]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
vs[3]=FullSimplify[matrixEigensystem[1][[2]][[3]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
vs[4]=FullSimplify[matrixEigensystem[1][[2]][[4]],Assumptions->{\[Chi]\[Element]Reals,\[Epsilon]\[Element]Reals,\[Omega]\[Element]Reals,\[Chi]>0,\[Epsilon]>0,\[Omega]>0}];
v1={v11,v12,v13,1};
v2={-v11,v12,-v13,1};
v3={v31,v32,v33,1};
v4={-v31,v32,-v33,1};
analSol[1]={c1[1]Exp[\[Lambda]1 t] v1,c1[2]Exp[-\[Lambda]1 t] v2,c1[3]Exp[\[Lambda]3 t] v3,c1[4]Exp[-\[Lambda]3 t] v4}\[Transpose];
analdf[1][t_]=(E^(-t (\[Lambda]1+\[Lambda]3)) (-E^(t (2 \[Lambda]1+\[Lambda]3)) (v11-v13) \[Lambda]3 c1[1]-E^(t \[Lambda]3) (v11-v13) \[Lambda]3 c1[2]-E^(t (\[Lambda]1+2 \[Lambda]3)) (v31-v33) \[Lambda]1 c1[3]-E^(t \[Lambda]1) (v31-v33) \[Lambda]1 c1[4]))/(\[Epsilon]^2 \[Lambda]1 \[Lambda]3)+c1[5];
analf[1][t_]=-((E^(-t (\[Lambda]1+\[Lambda]3)+t (2 \[Lambda]1+\[Lambda]3)) (v11-v13) c1[1])/(\[Epsilon]^2 \[Lambda]1^2))-(E^(-t (\[Lambda]1+\[Lambda]3)+t (\[Lambda]1+2 \[Lambda]3)) (v31-v33) c1[3])/(\[Epsilon]^2 \[Lambda]3^2)+(E^(-t (\[Lambda]1+\[Lambda]3)) ((E^(t \[Lambda]3) (v11-v13) \[Lambda]3 c1[2])/\[Lambda]1+(E^(t \[Lambda]1) (v31-v33) \[Lambda]1 c1[4])/\[Lambda]3))/(\[Epsilon]^2 \[Lambda]1 \[Lambda]3)+t c1[5]+c1[6];

analp[1][t_]=E^(t \[Lambda]1) v11 c1[1]-E^(-t \[Lambda]1) v11 c1[2]+E^(t \[Lambda]3) v31 c1[3]-E^(-t \[Lambda]3) v31 c1[4];
analdp[1][t_]=E^(t \[Lambda]1) v12 c1[1]+E^(-t \[Lambda]1) v12 c1[2]+E^(t \[Lambda]3) v32 c1[3]+E^(-t \[Lambda]3) v32 c1[4];
analm[1][t_]=E^(t \[Lambda]1) v13 c1[1]-E^(-t \[Lambda]1) v13 c1[2]+E^(t \[Lambda]3) v33 c1[3]-E^(-t \[Lambda]3) v33 c1[4];
analdm[1][t_]=E^(t \[Lambda]1) c1[1]+E^(-t \[Lambda]1) c1[2]+E^(t \[Lambda]3) c1[3]+E^(-t \[Lambda]3) c1[4];
analjp[1][t_]=E^(t \[Lambda]1) v12 c1[1]+E^(-t \[Lambda]1) v12 c1[2]+E^(t \[Lambda]3) v32 c1[3]+E^(-t \[Lambda]3) v32 c1[4]+(E^(-t (\[Lambda]1+\[Lambda]3)) (-E^(t (2 \[Lambda]1+\[Lambda]3)) (v11-v13) \[Lambda]3 c1[1]-E^(t \[Lambda]3) (v11-v13) \[Lambda]3 c1[2]-E^(t (\[Lambda]1+2 \[Lambda]3)) (v31-v33) \[Lambda]1 c1[3]-E^(t \[Lambda]1) (v31-v33) \[Lambda]1 c1[4]))/(\[Epsilon]^2 \[Lambda]1 \[Lambda]3)+c1[5];
anald\[Phi]d\[Phi][1][t_]=(1/(\[Epsilon]^2 \[Lambda]1 \[Lambda]3))(-E^(t \[Lambda]1) (v11-v13) \[Lambda]3 c1[1]+E^(-t \[Lambda]1) (-v11+v13) \[Lambda]3 c1[2]-E^(t \[Lambda]3) (v31-v33) \[Lambda]1 c1[3]+E^(-t \[Lambda]3) (-v31+v33) \[Lambda]1 c1[4]+\[Epsilon]^2 \[Lambda]1 \[Lambda]3 c1[5]) (-((E^(-t Conjugate[\[Lambda]1+\[Lambda]3]) Conjugate[E^(t (2 \[Lambda]1+\[Lambda]3)) (v11-v13) \[Lambda]3 c1[1]+E^(t \[Lambda]3) (v11-v13) \[Lambda]3 c1[2]+E^(t (\[Lambda]1+2 \[Lambda]3)) (v31-v33) \[Lambda]1 c1[3]+E^(t \[Lambda]1) (v31-v33) \[Lambda]1 c1[4]])/(\[Epsilon]^2 Conjugate[\[Lambda]1] Conjugate[\[Lambda]3]))+Conjugate[c1[5]]);
analjm[1][t_]=E^(t \[Lambda]1) c1[1]+E^(-t \[Lambda]1) c1[2]+E^(t \[Lambda]3) c1[3]+E^(-t \[Lambda]3) c1[4]-(E^(-t (\[Lambda]1+\[Lambda]3)) (-E^(t (2 \[Lambda]1+\[Lambda]3)) (v11-v13) \[Lambda]3 c1[1]-E^(t \[Lambda]3) (v11-v13) \[Lambda]3 c1[2]-E^(t (\[Lambda]1+2 \[Lambda]3)) (v31-v33) \[Lambda]1 c1[3]-E^(t \[Lambda]1) (v31-v33) \[Lambda]1 c1[4]))/(\[Epsilon]^2 \[Lambda]1 \[Lambda]3)-c1[5];
analConst[1][1]={{c1[1]->((1+E^(2 \[Lambda]3)) (1+v32) \[Epsilon]^2 \[Lambda]1^2 \[Lambda]3^2)/(2 (-2 E^\[Lambda]1 (1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]+4 E^(\[Lambda]1+\[Lambda]3) ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3])))),c1[2]->-((E^\[Lambda]1 (1+E^(2 \[Lambda]3)) (1+v32) \[Epsilon]^2 \[Lambda]1^2 \[Lambda]3^2)/(4 ((1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]-2 E^\[Lambda]3 ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3]))))),c1[3]->-(((1+E^(2 \[Lambda]1)) (1+v12) \[Epsilon]^2 \[Lambda]1^2 \[Lambda]3^2)/(2 (-2 E^\[Lambda]1 (1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]+4 E^(\[Lambda]1+\[Lambda]3) ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3]))))),c1[4]->-((E^(2 \[Lambda]3) (1+E^(2 \[Lambda]1)) (1+v12) \[Epsilon]^2 \[Lambda]1^2 \[Lambda]3^2)/(2 (-2 E^\[Lambda]1 (1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]+4 E^(\[Lambda]1+\[Lambda]3) ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3]))))),c1[5]->-(((1+E^(2 \[Lambda]1)) (1+E^(2 \[Lambda]3)) \[Lambda]1 \[Lambda]3 ((1+v12) (v31-v33) \[Lambda]1-(v11-v13) (1+v32) \[Lambda]3+(v12-v32) \[Epsilon]^2 \[Lambda]1 \[Lambda]3))/(2 (-2 E^\[Lambda]1 (1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]+4 E^(\[Lambda]1+\[Lambda]3) ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3]))))),c1[6]->((1+E^(2 \[Lambda]1)) (1+E^(2 \[Lambda]3)) \[Lambda]1 \[Lambda]3 ((1+v12) (v31-v33) \[Lambda]1-(v11-v13) (1+v32) \[Lambda]3+(v12-v32) \[Epsilon]^2 \[Lambda]1 \[Lambda]3))/(2 (-2 E^\[Lambda]1 (1+E^(2 \[Lambda]3)) \[Lambda]1 ((v11-v13) (1+v32)+(-v12+v32) \[Epsilon]^2 \[Lambda]1) \[Lambda]3^2 Cosh[\[Lambda]1]+4 E^(\[Lambda]1+\[Lambda]3) ((v11-v13) (1+v32) \[Lambda]3^2 Cosh[\[Lambda]3] Sinh[\[Lambda]1]+(1+v12) (v31-v33) \[Lambda]1^2 Cosh[\[Lambda]1] (\[Lambda]3 Cosh[\[Lambda]3]-Sinh[\[Lambda]3]))))}};

Begin[ "Private`"]
End[]
EndPackage[]
