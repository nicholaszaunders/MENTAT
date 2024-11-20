(* ::Package:: *)

BeginPackage["MENTAT`"];
Print["@MENTAT: Loading package..."]


isKet::usage="Returns True if argument is a ket object."
isKetState::usage="Return True if argument is a sum of ket objects."
isBra::usage="Returns True is argument is a bra object."
isBraState::usage="Returns True is argument is a sum of bra objects."
isDensity::usage="Returns True is argument is a density matrix object."
isDensityState::usage="Returns True is argument is a sum of density matrices."
isQuantumState::usage="Returns True if argument is a ket, sum of kets, bra, sum of bras, density matrix or sum of density matrices."
getPre::usage="Returns prefactor of the input ket, bra, or density matrix."
getNum::usage="Returns Fock basis sequence of the input ket or bra."
getNumDMLeft::usage="Returns Fock basis sequence of the ket half of input density matrix."
getNumDMRight::usage="Returns Fock basis sequence of the bra half of input density matrix."
makeCoherentState::usage="Takes complex amplitude \[Alpha] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the single-mode coherent state |\[Alpha]> up to cutoff."
makeTMSVSState::usage="Takes squeezing value \[Chi] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the two-mode squeezed vacuum state |\[Chi]> up to cutoff."
beamSplitter::usage="Implements beamsplitter functionality for a given ket-state x. Takes as input a single ket x, the beamsplitter transmissivity \[Tau], and the indices of the two modes to be mixed. Returns the ket state x' after the unitary beamsplitter is applied."
partialTrace::usage="Partial trace calculator. Takes as input a density matrix x of a state with n modes and a list tracedModesList of length < n describing the indices of the modes to be traced out. Returns a density matrix x' corresponding to x with the specified modes traced out."
generateNDArray::usage="Generates a list containing the basis Fock states of a Hilbert space corresponding to a system with number of modes numModes and maximum Fock-state argument cutoff. The list is of length (cutoff+1)^numModes."
generateDMMatrixRepresentation::usage="Takes a given density matrix x and the Hilbert space characteristics numModes, cutoff and produces a numeric matrix isomorphic to that density matrix in the basis space {|n1,n2,...,nNumModes>} for \!\(\*FormBox[\(\*SubscriptBox[\(n\), \(i\)] < cutoff\),
TraditionalForm]\). Returns square matrix of dimension (cutoff+1)^numModes."
calculateNumericalMatrixFidelity::usage="Calculate the quantum fidelity of two positive semidefinite matrices x, y assuming they correspond to normalised and valid quantum states."
vonNeumannEntropy::usage="Calculate the von Neumann entropy S(\[Rho]) of a density matrix state \[Rho]."


Begin["`Private`"];


(* ::Section:: *)
(*Initialisation*)


(*Identifier function for the 'ket' object. Returns True if the object is a ket.
 The algebra uses Mathematica's built-in Times to deal with scalar multiplication, so I define a ket as any sequence with the built-in header Ket or any sequence with the header Ket multiplied by some arbitrary quantity.*)
 isKet[expr_]:=
MatchQ[
expr,
__Ket|Times[_,__Ket]
]


(*Identifier function for the 'ket state' object. Returns True if the object is a sum of kets.
I don't know why I need HoldPattern here, but I definitely do.*)
isKetState[expr_]:=
MatchQ[
expr,
HoldPattern[\!\(\*
TagBox[
StyleBox[
RowBox[{"Plus", "[", 
RowBox[{"__", "?", "isKet"}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)]
]


(*Identifier function for the 'bra' object. Returns True if the object is a bra.
Uses the same rules as the ket object but for the built-in header Bra.*)
isBra[expr_]:=
MatchQ[
expr,
__Bra|Times[_,__Bra]
]


(*See isKetState.*)
isBraState[expr_]:=
MatchQ[
expr,
HoldPattern[\!\(\*
TagBox[
StyleBox[
RowBox[{"Plus", "[", 
RowBox[{"__", "?", "isBra"}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)]
]


(*Identifier function for the 'density matrix' object. Returns true if expr is a ket object vector 'small-circle'-multiplied by a bra object.
We use the SmallCircle operator \[SmallCircle] / \[SmallCircle] to fill the definition of a density state class.
Annoyingly, even though technically we should be able to use our vector multiplication operator \[CenterDot] to fill both the inner and outer product, this doesn't play well with Mathematica because using the same operator for both operations and for 'definitions' of the density matrix causes recursion. I could use the same operator if I found a way to manually halt recursion when evaluating, but using two different operators is easiest for now.*)
isDensity[expr_]:=
MatchQ[
expr,
SmallCircle[__Ket,__Bra]|Times[_,SmallCircle[__Ket,__Bra]]
]


(*See isKetState.*)
isDensityState[expr_]:=
MatchQ[
expr,
HoldPattern[\!\(\*
TagBox[
StyleBox[
RowBox[{"Plus", "[", 
RowBox[{"__", "?", "isDensity"}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)]
]


(*Identifier function for any object that is a valid quantum state,i . e . any expression that is either a ket,bra,density matrix or sum of the respective objects .*)
isQuantumState[expr_]:=
MatchQ[
expr,
(_?isKet)|(_?isBra)|(_?isDensity)|(_?isKetState)|(_?isBraState)|(_?isDensityState)
]


(*Returns the prefactor pre of expr when expr is a ket object,bra object,or density matrix object . Different cases have to be put in to account for whether the object has a prefactor or not . We do this because using Mathematica's built-in Times to deal with scalar multiplication means objects multiplied by 1 or 0 will have a different code structure than those multiplied by an arbitrary quantity .*)
getPre[expr_?isKet|expr_?isBra|expr_?isDensity]:=
Which[
MatchQ[expr,__Ket],
	1,
MatchQ[expr,Times[_,__Ket]],
	expr/.Times[pre_,Ket[__]]->pre,
MatchQ[expr,__Bra],
	1,
MatchQ[expr,Times[_,__Bra]],
	expr/.Times[pre_,Bra[__]]->pre,
MatchQ[expr,SmallCircle[Ket[__],Bra[__]]],
	1,
MatchQ[expr,Times[_,SmallCircle[__Ket,__Bra]]],
	expr/.Times[pre_,SmallCircle[Ket[__],Bra[__]]]->pre
]


(*Returns the Sequence num, describing the argument or mode numbers, of expr when expr is a ket or bra.
See getPre.*)
getNum[expr_?isKet|expr_?isBra]:=
Which[
MatchQ[expr,__Ket],
	expr/.Ket[num__]->Sequence[num],
MatchQ[expr,Times[_,__Ket]],
	expr/.Times[_,Ket[num__]]->num,
MatchQ[expr,__Bra],
	expr/.Bra[num__]->Sequence[num],
MatchQ[expr,Times[_,__Bra]],
	expr/.Times[_,Bra[num__]]->num
]


(*Retrieves Sequence describing the mode numbers num of the ket state making up a density matrix object expr.*)
getNumDMLeft[expr_?isDensity]:=
Which[
MatchQ[expr,SmallCircle[__Ket,__Bra]],
	expr/.SmallCircle[Ket[num__],_]->num,
MatchQ[expr,Times[_,SmallCircle[__Ket,__Bra]]],
	expr/.Times[_,SmallCircle[Ket[num__],_]]->num
]


(*Retrieves Sequence describing the mode numbers num of the bra state making up a density matrix object expr.*)
getNumDMRight[expr_?isDensity]:=
Which[
MatchQ[expr,SmallCircle[__Ket,__Bra]],
	expr/.SmallCircle[_,Bra[num__]]->num,
MatchQ[expr,Times[_,SmallCircle[__Ket,__Bra]]],
	expr/.Times[_,SmallCircle[_,Bra[num__]]]->num
]


(* ::Section:: *)
(*Inner product \[CenterDot]*)


(*Define operation of vector multiplication CenterDot for the vector inner product of two pure states. Returns scalar.
If x is a bra and y a ket, then return 0 if num of x != num of y, and return the product of the prefactors otherwise.*)
CenterDot[x_?isBra,y_?isKet]:=
If[
List[getNum[x]]===List[getNum[y]],
getPre[x]*getPre[y],
0
]


(*Define vector multiplication for vector outer product of two pure states. Returns density matrix.
If x is a ket and y a bra, then 
	If both prefactors are 1 (recalling that we define density matrix objects explicitly via the CenterDot of a ket and bra) do nothing
	If a prefactor is non-unity, simplify by bringing it to the front.
Note that we have to use a more general pattern-matching strategy (i.e. not isBra or isKet) because we want to avoid having to explicitly define the edge case of a pure bra times a pure ket; this causes infinite recursion.*)
SmallCircle[Times[preKet_,Ket[numKet__]],Bra[numBra__]]:=
Times[preKet,SmallCircle[Ket[numKet],Bra[numBra]]]

SmallCircle[Ket[numKet__],Times[preBra_,Bra[numBra__]]]:=
Times[preBra,SmallCircle[Ket[numKet],Bra[numBra]]]

SmallCircle[Times[preKet_,Ket[numKet__]],Times[preBra_,Bra[numBra__]]]:=
Times[preKet,preBra,SmallCircle[Ket[numKet],Bra[numBra]]]


(*Define vector multiplication for inner product between two density matrices. Returns density matrix.*)
CenterDot[x_?isDensity,y_?isDensity]:=
If[
List[getNumDMRight[x]]===List[getNumDMLeft[y]],
getPre[x]*getPre[y]*SmallCircle[Ket[getNumDMLeft[x]],Bra[getNumDMRight[y]]],
0
]


(*Define vector multiplication for inner product between left-bra and right-DM. Returns bra.*)
CenterDot[x_?isBra,y_?isDensity]:=
If[
List[getNum[x]]===List[getNumDMLeft[y]],
getPre[x]*getPre[y]*Bra[getNumDMRight[y]],
0
]


(*Define vector multiplication for inner product between left-DM and right-ket. Returns ket.*)
CenterDot[x_?isDensity,y_?isKet]:=
If[
List[getNumDMRight[x]]===List[getNum[y]],
getPre[x]*getPre[y]*Ket[getNumDMLeft[x]],
0
]


(*Extend inner product to operations with scalars.
Cover distributivity with scalars.*)
CenterDot[x_?isKet|x_?isBra|x_?isDensity,Except[y_?isQuantumState]]:=
y*x

CenterDot[Except[x_?isQuantumState],y_?isKet|y_?isBra|y_?isDensity]:=
x*y

CenterDot[x_?isKetState|x_?isBraState|x_?isDensityState,Except[y_?isQuantumState]]:=
Sum[
y*x[[i]],
{i,Length[x]}
]

CenterDot[Except[x_?isQuantumState],y_?isKetState|y_?isBraState|y_?isDensityState]:=
Sum[
x*y[[i]],
{i,Length[y]}
]


(*Extend inner product to bra and ket states (distributivity).
Additional logic needs to be put in for the degenerate case where a single ket or bra is multplied by a ket state or bra state. This probably isn't very efficient, since the amount of different cases grows exponentially. Fortunately we only really need to consider bras, kets, and density matrices. Another important reason for doing it this way is that we can explicitly exclude unwanted scenarios, like ket\[CenterDot]DM or DM\[CenterDot]bra.*)
CenterDot[x_?isBraState,y_?isKetState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isBraState,y_?isKet]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]
CenterDot[x_?isBra,y_?isKetState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]

CenterDot[x_?isKetState,y_?isBraState]:=
Sum[
SmallCircle[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isKet,y_?isBraState]:=
Sum[
SmallCircle[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isKetState,y_?isBra]:=
Sum[
SmallCircle[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isDensityState,y_?isBraState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isDensity,y_?isBraState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isDensityState,y_?isBra]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isBraState,y_?isDensityState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isBraState,y_?isDensity]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]
CenterDot[x_?isBra,y_?isDensityState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]

CenterDot[x_?isKetState,y_?isDensityState]:=
Sum[
SmallCircle[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isKet,y_?isDensityState]:=
Sum[
SmallCircle[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isKetState,y_?isDensity]:=
Sum[
SmallCircle[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isDensityState,y_?isKetState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isDensityState,y_?isKet]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]
CenterDot[x_?isDensity,y_?isKetState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]

CenterDot[x_?isDensityState,y_?isDensityState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isDensityState,y_?isDensity]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]
CenterDot[x_?isDensity,y_?isDensityState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]


(*Extend inner product to expectation-value type products, where the multiplication is of the form \[LeftAngleBracket]a|\[CenterDot]b\[CenterDot]|c\[RightAngleBracket].*)
CenterDot[x_?isBra|x_?isBraState,z_,y_?isKet|y_?isKetState]:=
CenterDot[x,CenterDot[z,y]]


(* ::Section:: *)
(*Outer product \[CircleTimes]*)


(*Base definition of the tensor product. 
We restrict operation to kets since I can't think of any reason why bras or density matrices would need to be directly tensor operated on.*)
CircleTimes[x_?isKet,y_?isKet]:=
getPre[x]*getPre[y]Ket[getNum[x],getNum[y]]


(*Sets associativity, i.e. Ket[1]\[CircleTimes]Ket[1]\[CircleTimes]Ket[1] = Ket[1,1,1].
For a given n-tensor product, the product is replaced by a list, and in that list each Ket object is replaced by the Sequence of Fock-basis arguments. The list is then transformed
back into a Ket object.*)

CircleTimes[x__?isKet]:=
Times@@(List[x]/.y_?isKet->getPre[y])*Ket@@(List[x]/.y_?isKet->getNum[y])


(*Defining compatibility with scalar (0-dimension) quantities*)
CircleTimes[x_?isKet,Except[y_?isKet|y_?isKetState]]:=
getPre[x]*y Ket[getNum[x]]

CircleTimes[Except[x_?isKet|x_?isKetState],y_?isKet]:=
getPre[y]*x Ket[getNum[y]]

CircleTimes[Except[x_?isKet|x_?isKetState],Except[y_?isKet|y_?isKetState]]:=
x*y


(*Defining distributivity in the same manner as above.*)
CircleTimes[x__?isKetState,y__?isKetState]:=
Sum[
CircleTimes[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CircleTimes[x__?isKet,y__?isKetState]:=
Sum[
CircleTimes[x,y[[j]]],
{j,Length[y]}
]
CircleTimes[x__?isKetState,y__?isKet]:=
Sum[
CircleTimes[x[[i]],y],
{i,Length[x]}
]


(* ::Section:: *)
(*Functionality*)


(*Returns a ket state describing the coherent state Ket[\[Alpha]]for a complex amplitude \[Alpha] and maximum Fock-state element cutoff.*)
makeCoherentState[\[Alpha]_,cutoff_]:=
Sum[E^(-(Abs[\[Alpha]]^2/2)) \[Alpha]^n/\[Sqrt](n!) Ket[n],{n,0,cutoff}]


(*Returns a two-mode squeezed vacuum state Ket[\[Chi]] for squeezing \[Chi] = tanh(r) and maximum Fock-state element cutoff.*)
makeTMSVSState[\[Chi]_,cutoff_]:=
Sum[\[Sqrt](1-\[Chi]^2) \[Chi]^n Ket[n,n],{n,0,cutoff}]


(*Implements beamsplitter functionality for a given ket-state x. Takes as input a single ket x, the beamsplitter transmissivity \[Tau], and the indices of the two modes to be mixed. Returns the ket state x' after the unitary beamsplitter is applied.*)
beamSplitter[x_?isKet,\[Tau]_,mode1_,mode2_]:=
Module[{modeList,n,m},
modeList=List[getNum[x]];
n=List[getNum[x]][[mode1]];
m=List[getNum[x]][[mode2]];
Expand[
Sum[
getPre[x](-1)^i (Sqrt[(i+j)!] Sqrt[m!] Sqrt[n!] Sqrt[(n-i+m-j)!])/(i! j! (m-j)! (n-i)!) \[Sqrt](\[Tau])^(n-i+j) \[Sqrt](1-\[Tau])^(i+m-j)ReplaceAt[Ket[getNum[x]],{modeList[[mode1]]->n-i+m-j,modeList[[mode2]]->i+j},{{mode1},{mode2}}],
{i,0,n},
{j,0,m}
]
]
]


(*Defines distributivity for the beamsplitter operation.*)
beamSplitter[x__?isKetState,\[Tau]_,mode1_,mode2_]:=
Sum[
beamSplitter[x[[i]],\[Tau],mode1,mode2],
{i,Length[x]}
]


(*Defines operation of the conjugate tranpose operation ^Dagger[*]. For an input quantum state x, conjugate the prefactor and transpose the quantum state appropriately.*) 
SuperDagger[x_?isKet]:=
Conjugate[getPre[x]]*Bra[getNum[x]]

SuperDagger[x_?isBra]:=
Conjugate[getPre[x]]*Ket[getNum[x]]

SuperDagger[x_?isDensity]:=
Conjugate[getPre[x]]*SmallCircle[Ket[getNumDMRight[x]],Bra[getNumDMLeft[x]]]


(*Define distributivity of the conjugate transpose.*)
SuperDagger[x_?isKetState]:=
Sum[
SuperDagger[x[[i]]],
{i,Length[x]}
]

SuperDagger[x_?isBraState]:=
Sum[
SuperDagger[x[[i]]],
{i,Length[x]}
]

SuperDagger[x_?isDensityState]:=
Sum[
SuperDagger[x[[i]]],
{i,Length[x]}
]


(*Partial trace calculator. Takes as input a density matrix x of a state with n modes and a list tracedModesList of length < n describing the indices of the modes to be traced out. Returns a density matrix x' corresponding to x with the specified modes traced out.
Note: this is less of a mathematically rigorous definition and more of a computational one. Faster and easier, but it's a shame to lose the mathematical backing.*)
partialTrace[x_?isDensity,tracedModesList_]:=
Module[
{numLeft, numRight,tracedNumLeft,tracedNumRight},
numLeft=List[getNumDMLeft[x]];
numRight=List[getNumDMRight[x]];
tracedNumLeft=Delete[numLeft,List/@tracedModesList];
tracedNumRight=Delete[numRight,List/@tracedModesList];
If[
numLeft[[tracedModesList]]===numRight[[tracedModesList]],
getPre[x]*SmallCircle[Ket[tracedNumLeft/.List:>Sequence],Bra[tracedNumRight/.List:>Sequence]],
0
]
]


(*Define distributivity of the partial trace.*)
partialTrace[x_?isDensityState,tracedModesList_]:=
Sum[
partialTrace[x[[i]],tracedModesList],
{i,Length[x]}
]


(*Generates a list containing the basis Fock states of a Hilbert space corresponding to a system with number of modes numModes and maximum Fock-state argument cutoff.
The list is of length (cutoff+1)^numModes.*)
generateNDArray[numModes_,cutoff_]:=
Module[
{basisVectorSum},
basisVectorSum=Sum[Ket[i],{i,0,cutoff}];
Do[
basisVectorSum=(basisVectorSum)\[CircleTimes](Sum[Ket[i],{i,0,cutoff}]),
{j,1,numModes-1}
];
Return[basisVectorSum/.Plus->List]
]


(*Takes a given density matrix x and the Hilbert space characteristics numModes, cutoff and produces a numeric matrix isomorphic to that density matrix in the basis space {|n1,n2,...,nNumModes>} for Subscript[n, i]<cutoff. Returns square matrix of dimension (cutoff+1)^numModes.*)
generateDMMatrixRepresentation[x_?isDensityState,numModes_,cutoff_]:=
Module[
{rhoMatrix,basisList},
rhoMatrix=ConstantArray[0,{(cutoff+1)^numModes,(cutoff+1)^numModes}];
basisList=generateNDArray[numModes,cutoff];
Do[
rhoMatrix[[i,j]]=SuperDagger[(basisList[[i]])]\[CenterDot](x\[CenterDot]basisList[[j]]),
{i,1,(cutoff+1)^numModes},
{j,1,(cutoff+1)^numModes}
];
Return[rhoMatrix]
]


(*Calculate the quantum fidelity of two positive semidefinite matrices x, y assuming they correspond to normalised and valid quantum states.*)
calculateNumericalMatrixFidelity[x_,y_]:=
Tr[MatrixPower[MatrixPower[x,1/2] . y . MatrixPower[x,1/2],1/2]]^2


(*Calculate the von Neumann entropy of a given density matrix state.*)
vonNeumannEntropy[x_?isDensityState,numModes_,cutoff_]:=
Module[
{numericalMatrix,eigenvalues},
numericalMatrix=generateDMMatrixRepresentation[x,numModes,cutoff];
eigenvalues=Eigenvalues[numericalMatrix];
Sum[
If[
eigenvalues[[i]]==0,
0,
-eigenvalues[[i]]*Log2[eigenvalues[[i]]]
],
{i,1,Length[eigenvalues]}
]
]


Print["@MENTAT: All functions loaded. 
Welcome to MENTAT, a CAS system for Fock-state calculations in bosonic optical systems! Copyright N. Zaunders, University of Queensland, 2024. All rights reserved."]


End[];


EndPackage[];
