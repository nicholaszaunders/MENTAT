(* ::Package:: *)

(* ::Text:: *)
(*MENTAT package + documentation*)
(*Nov 2024*)
(*N. Zaunders, University of Queensland School of Mathematics and Physics*)
(*To do:*)
(*- Might be nice to do some functionality where you can specify a given norm P, i.e. 99.999, and the machine chooses the cutoff automatically such that the elements contained within the truncated state make up >= P of the continuous state.*)
(*- Change associativity to operate from right to left instead of left to right for better integration with QM?*)
(*- Now that I've solved the operator formalism, maybe switch over some of the functionality / definitions to use operators instead of hardcoding?*)


BeginPackage["MENTAT`"];
Print["@MENTAT: Clearing kernel..."]
Quit
Print["@MENTAT: Loading package..."]


isKet::usage="Returns True if argument is a ket object."
isKetState::usage="Return True if argument is a sum of ket objects."
isBra::usage="Returns True is argument is a bra object."
isBraState::usage="Returns True is argument is a sum of bra objects."
isDensity::usage="Returns True is argument is a density matrix object."
isDensityState::usage="Returns True is argument is a sum of density matrices."
isQuantum::usage="Returns True if argument is a ket, sum of kets, bra, sum of bras, density matrix or sum of density matrices."
getPre::usage="Returns prefactor of the input ket, bra, or density matrix."
getNum::usage="Returns Fock basis sequence of the input ket or bra."
getNumDMLeft::usage="Returns Fock basis sequence of the ket half of input density matrix."
getNumDMRight::usage="Returns Fock basis sequence of the bra half of input density matrix."

isIdentity::usage="Returns True if argument is the identity operator \[ScriptCapitalI]."
isCreationOperator::usage="Returns True if argument is a creation operator of the form \!\(\*SuperscriptBox[SubscriptBox[OverscriptBox[\(a\), \(^\)], \(i\)], \(\[Dagger]\)]\). Can have prefactors and be raised to any power."
isAnnihilationOperator::usage="Returns True if argument is an annihilation operator of the form \!\(\*SubscriptBox[OverscriptBox[\(a\), \(^\)], \(i\)]\). Can have prefactors and be raised to any power."
isOperatorComposition::usage="Returns True is argument is a CenterDot composition of creation or annihilation operators."
isOperator::usage="Returns True if argument is a creation or annihilation operator or composition thereof."
isOperatorSum::usage="Returns True if argument is a sum of creation or annihilation operator raised to any power."
getCreateAnnihilateOperatorMode::usage="Returns mode index 'mode' of the input creation or annihilation operator \!\(\*SubscriptBox[OverscriptBox[\(a\), \(^\)], \(i\)]\)."
getCreateAnnihilateOperatorPower::usage="Returns power pow of the input creation or annihilation operator."
getCreateAnnihilateOperatorPre::usage="Returns prefactor pre of the input creation or annihilation operator."

makeCoherentState::usage="Takes complex amplitude \[Alpha] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the single-mode coherent state |\[Alpha]> up to cutoff."
makeTMSVSState::usage="Takes squeezing value \[Chi] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the two-mode squeezed vacuum state |\[Chi]> up to cutoff."
beamSplitter::usage="Implements beamsplitter functionality for a given ket-state x. Takes as input a single ket x, the beamsplitter transmissivity \[Tau], and the indices of the two modes to be mixed. Returns the ket state x' after the unitary beamsplitter is applied."
partialTrace::usage="Partial trace calculator. Takes as input a density matrix x of a state with n modes and a list tracedModesList of length < n describing the indices of the modes to be traced out. Returns a density matrix x' corresponding to x with the specified modes traced out."
generateNDArray::usage="Generates a list containing the basis Fock states of a Hilbert space corresponding to a system with number of modes numModes and maximum Fock-state argument cutoff. The list is of length (cutoff+1)^numModes."
generateDMMatrixRepresentation::usage="Takes a given density matrix x and the Hilbert space characteristics numModes, cutoff and produces a numeric matrix isomorphic to that density matrix in the basis space {|n1,n2,...,nNumModes>} for \!\(\*FormBox[\(\*SubscriptBox[\(n\), \(i\)] < cutoff\),
TraditionalForm]\). Returns square matrix of dimension (cutoff+1)^numModes."
calculateNumericalMatrixFidelity::usage="Calculate the quantum fidelity of two positive semidefinite matrices x, y assuming they correspond to normalised and valid quantum states."
vonNeumannEntropy::usage="Calculate the von Neumann entropy S(\[Rho]) of a density matrix state \[Rho]."
fullTrace::usage="Calculates the trace of an input density matrix."


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


(*Identifier function for any object that is a valid quantum object, i.e. any expression that is either a ket, bra, density matrix, sum of the respective objects, operator, etc.
Essentially acts as a reverse filter to identify scalars, which can be numerical, symbolic, etc.*)
isQuantum[expr_]:=
MatchQ[
expr,
_?isKet|_?isBra|_?isDensity|_?isKetState|_?isBraState|_?isDensityState|_?isOperator|_?isOperatorSum
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
CenterDot[x_?isKet|x_?isBra|x_?isDensity,Except[y_?isQuantum]]:=
y*x

CenterDot[Except[x_?isQuantum],y_?isKet|y_?isBra|y_?isDensity]:=
x*y

CenterDot[x_?isKetState|x_?isBraState|x_?isDensityState,Except[y_?isQuantum]]:=
Sum[
y*x[[i]],
{i,Length[x]}
]

CenterDot[Except[x_?isQuantum],y_?isKetState|y_?isBraState|y_?isDensityState]:=
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


(*Extend inner product to be associative for expectation-value type products, where the multiplication is of the form \[LeftAngleBracket]a|\[CenterDot]b\[CenterDot]|c\[RightAngleBracket].*)
CenterDot[x_?isBra|x_?isBraState,z_,y_?isKet|y_?isKetState]:=
CenterDot[x,CenterDot[z,y]]


(* ::Section:: *)
(*Outer product \[CircleTimes]*)


(*Base definition of the tensor product. *)
CircleTimes[x_?isKet,y_?isKet]:=
getPre[x]*getPre[y]*Ket[getNum[x],getNum[y]]

CircleTimes[x_?isBra,y_?isBra]:=
getPre[x]*getPre[y]*Bra[getNum[x],getNum[y]]

CircleTimes[x_?isDensity,y_?isDensity]:=
getPre[x]*getPre[y]*SmallCircle[Ket[getNumDMLeft[x],getNumDMLeft[y]],Bra[getNumDMRight[x],getNumDMRight[y]]]


(*Defining compatibility with scalar (0-dimension) quantities.*)
CircleTimes[x_?isKet|x_?isBra|x_?isDensity,Except[y_?isQuantum]]:=
y*x

CircleTimes[Except[x_?isQuantum],y_?isKet|y_?isBra|y_?isDensity]:=
x*y

CircleTimes[x_?isKetState|x_?isBraState|x_?isDensityState,Except[y_?isQuantum]]:=
Sum[
y*x[[i]],
{i,Length[x]}
]

CircleTimes[Except[x_?isQuantum],y_?isKetState|y_?isBraState|y_?isDensityState]:=
Sum[
x*y[[i]],
{i,Length[y]}
]


(*Sets associativity, i.e. Ket[1]\[CircleTimes]Ket[1]\[CircleTimes]Ket[1] = Ket[1,1,1].*)

CircleTimes[x__/;Length[List[x]]>=3]:=
Module[
	{args,len,iterProduct},
	args=List[x];
	len=Length[args];
Do[
	iterProduct=CircleTimes[args[[1]],args[[2]]];
	args=Drop[args,1];
	args[[1]]=iterProduct,
	{i,len-2}
];
Return[CircleTimes[args[[1]],args[[2]]]]
]


(*Defining distributivity in the same manner as above.*)
CircleTimes[x_?isKetState,y_?isKetState]:=
Sum[
CircleTimes[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CircleTimes[x_?isKetState,y_?isKet]:=
Sum[
CircleTimes[x[[i]],y],
{i,Length[x]}
]
CircleTimes[x_?isKet,y_?isKetState]:=
Sum[
CircleTimes[x,y[[j]]],
{j,Length[y]}
]

CircleTimes[x_?isBraState,y_?isBraState]:=
Sum[
CircleTimes[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CircleTimes[x_?isBraState,y_?isBra]:=
Sum[
CircleTimes[x[[i]],y],
{i,Length[x]}
]
CircleTimes[x_?isBra,y_?isBraState]:=
Sum[
CircleTimes[x,y[[j]]],
{j,Length[y]}
]

CircleTimes[x_?isDensityState,y_?isDensityState]:=
Sum[
CircleTimes[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CircleTimes[x_?isDensityState,y_?isDensity]:=
Sum[
CircleTimes[x[[i]],y],
{i,Length[x]}
]
CircleTimes[x_?isDensity,y_?isDensityState]:=
Sum[
CircleTimes[x,y[[j]]],
{j,Length[y]}
]


(* ::Section:: *)
(*Operators*)


(* ::Subsection:: *)
(*Identity operator*)


(*Defining matching function for identity operator on a mode.*)
(* !!! UNFINISHED !!! *)
(*For some reason, the below function works when defined in a notebook but doesn't work when defined here. Revisit later.*)
(*isIdentity[expr_]:=
SameQ[expr,\[DoubleStruckCapitalI]]*)


(*Defining vector multiplication dot product for identity operator.*)
(* !!! UNFINISHED !!! *)
(*CenterDot[x_?isQuantum,y_?isIdentity]:=
x

CenterDot[x_?isIdentity,y_?isQuantum]:=
y*)


(* ::Subsection:: *)
(*Creation and annihilation operators*)


(*Define identifier function for creation operators. Returns True if expr is a creation operator.*)
isCreationOperator[expr_]:=
MatchQ[
expr,
	 SuperDagger[Subscript[OverHat[a_Symbol],_Integer]]
	|Times[_,SuperDagger[Subscript[OverHat[a_Symbol],_Integer]]]
	|Power[SuperDagger[Subscript[OverHat[a_Symbol],_Integer]],_Integer]
	|Times[_,Power[SuperDagger[Subscript[OverHat[a_Symbol],_Integer]],_Integer]]
]

(*Define identifier function for powers of annihilation operators. Returns True if expr is an annihilation operator.*)
isAnnihilationOperator[expr_]:=
MatchQ[
expr,
	 Subscript[OverHat[a_Symbol],_Integer]
	|Times[_,Subscript[OverHat[a_Symbol],_Integer]]
	|Power[Subscript[OverHat[a_Symbol],_Integer],_Integer]
	|Times[_,Power[Subscript[OverHat[a_Symbol],_Integer],_Integer]]
]

(*Define identifier function for a composition of Fock state operators, i.e. Subscript[Overscript[a, ^], i]^\[Dagger]\[CenterDot]Subscript[Overscript[a, ^], j]^\[Dagger]\[CenterDot]..., etc. Returns True if expr is an arbitrary series of creation or annihilation operators
CenterDot'd with each other.
Works by checking truth of isCreation/isAnnihilation for each element in the Sequence passed to CenterDot. The Sequence of booleans is then checked via And; if every element returns true, 
then it must be a composition of operators.*)
isOperatorComposition[expr_]:=
MatchQ[
expr,
CenterDot[x__/;(And@@(MatchQ[_?isCreationOperator|_?isAnnihilationOperator]/@(List[x])))]
]

(*Define identifier function for any type of Fock state operator. Returns True if expr is a creation or annihilation operator or composition thereof.*)
isOperator[expr_]:=
MatchQ[
expr,
_?isCreationOperator|_?isAnnihilationOperator|_?isOperatorComposition
]

(*Define identifier function for any sum of Fock state operators. Returns True if expr is a sum of creation or annihilation operators or composition thereof.*)
isOperatorSum[expr_]:=
MatchQ[
expr,
HoldPattern[\!\(\*
TagBox[
StyleBox[
RowBox[{"Plus", "[", 
RowBox[{"__", "?", "isOperator"}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)]
]

(*Returns the mode index that the creation or annihilation operator x should operate on.*)
getCreateAnnihilateOperatorMode[expr_?isCreationOperator|expr_?isAnnihilationOperator]:=
Which[
isCreationOperator[expr],
	Which[
	MatchQ[expr,SuperDagger[Subscript[OverHat[_],mode_Integer]]],
		expr/.SuperDagger[Subscript[OverHat[_],mode_Integer]]->mode,
		
	MatchQ[expr,Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]],
		expr/.Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]->mode,
		
	MatchQ[expr,Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]],
		expr/.Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]->mode,
		
	MatchQ[expr,Times[_,Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]],
		expr/.Times[_,Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]->mode
	],
isAnnihilationOperator[expr],
	Which[
	MatchQ[expr,Subscript[OverHat[_],mode_Integer]],
		expr/.Subscript[OverHat[_],mode_Integer]->mode,
		
	MatchQ[expr,Times[_,Subscript[OverHat[_],mode_Integer]]],
		expr/.Times[_,Subscript[OverHat[_],mode_Integer]]->mode,
		
	MatchQ[expr,Power[Subscript[OverHat[_],mode_Integer],_]],
		expr/.Power[Subscript[OverHat[_],mode_Integer],_]->mode,
		
	MatchQ[expr,Times[_,Power[Subscript[OverHat[_],mode_Integer],_]]],
		expr/.Times[_,Power[Subscript[OverHat[_],mode_Integer],_]]->mode
	]
]

(*Returns the power that the creation or annihilation operator expr is raised to.*)
getCreateAnnihilateOperatorPower[expr_?isCreationOperator|expr_?isAnnihilationOperator]:=
Which[
isCreationOperator[expr],
	Which[
	MatchQ[expr,SuperDagger[Subscript[OverHat[_],mode_Integer]]],
		1,
	MatchQ[expr,Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]],
		1,
	MatchQ[expr,Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]],
		expr/.Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]->pow,
	MatchQ[expr,Times[_,Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]],
		expr/.Times[_,Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]->pow
	],
isAnnihilationOperator[expr],
	Which[
	MatchQ[expr,Subscript[OverHat[_],_]],
		1,
	MatchQ[expr,Times[_,Subscript[OverHat[_],_]]],
		1,
	MatchQ[expr,Power[Subscript[OverHat[_],_],pow_Integer]],
		expr/.Power[Subscript[OverHat[_],_],pow_Integer]->pow,
	MatchQ[expr,Times[_,Power[Subscript[OverHat[_],_],pow_Integer]]],
		expr/.Times[_,Power[Subscript[OverHat[_],_],pow_Integer]]->pow
	]
]

(*Returns the prefactor of the creation or annihilation operator.*)
getCreateAnnihilateOperatorPre[expr_?isCreationOperator|expr_?isAnnihilationOperator]:=
Which[
isCreationOperator[expr],
	Which[
	MatchQ[expr,SuperDagger[Subscript[OverHat[_],_]]],
		1,
	MatchQ[expr,Times[pre_,SuperDagger[Subscript[OverHat[_],_]]]],
		expr/.Times[pre_,SuperDagger[Subscript[OverHat[_],_]]]->pre,
	MatchQ[expr,Power[SuperDagger[Subscript[OverHat[_],_]],_]],
		1,
	MatchQ[expr,Times[pre_,Power[SuperDagger[Subscript[OverHat[_],_]],_]]],
		expr/.Times[pre_,Power[SuperDagger[Subscript[OverHat[_],_]],_]]->pre
	],
isAnnihilationOperator[expr],
	Which[
	MatchQ[expr,Subscript[OverHat[_],_]],
		1,
	MatchQ[expr,Times[pre_,Subscript[OverHat[_],_]]],
		expr/.Times[pre_,Subscript[OverHat[_],_]]->pre,
	MatchQ[expr,Power[Subscript[OverHat[_],_],_]],
		1,
	MatchQ[expr,Times[pre_,Power[Subscript[OverHat[_],_],_]]],
		expr/.Times[pre_,Power[Subscript[OverHat[_],_],_]]->pre
	]
]


(*Define raising and lowering operations. We re-use the operator \[CenterDot] to signify interaction with an operator.
To give multimodal functionality, each operator is equipped with an integer subscript which corresponds to the mode it should operate on.*)
CenterDot[x_?isCreationOperator,y_?isKet]:=
Module[
	{operatorMode, modeList},
	modeList=List[getNum[y]];
	operatorMode=getCreateAnnihilateOperatorMode[x];
Sqrt[modeList[[operatorMode]]+1]*getPre[y]*getCreateAnnihilateOperatorPre[x]*ReplaceAt[Ket[getNum[y]],{modeList[[operatorMode]]->modeList[[operatorMode]]+1},{operatorMode}]
]

CenterDot[x_?isAnnihilationOperator,y_?isKet]:=
Module[
	{operatorMode, modeList},
	modeList=List[getNum[y]];
	operatorMode=getCreateAnnihilateOperatorMode[x];
If[
	modeList[[operatorMode]]=!=0,
	Sqrt[modeList[[operatorMode]]]*getPre[y]*getCreateAnnihilateOperatorPre[x]*ReplaceAt[Ket[getNum[y]],{modeList[[operatorMode]]->modeList[[operatorMode]]-1},{operatorMode}],
	0
]
]

CenterDot[x_?isBra,y_?isCreationOperator]:=
Module[
	{operatorMode, modeList},
	modeList=List[getNum[x]];
	operatorMode=getCreateAnnihilateOperatorMode[y];
If[
	modeList[[operatorMode]]=!=0,
	Sqrt[modeList[[operatorMode]]]*getPre[x]*getCreateAnnihilateOperatorPre[y]*ReplaceAt[Bra[getNum[x]],{modeList[[operatorMode]]->modeList[[operatorMode]]-1},{operatorMode}],
	0
]
]

CenterDot[x_?isBra,y_?isAnnihilationOperator]:=
Module[
	{operatorMode, modeList},
	modeList=List[getNum[x]];
	operatorMode=getCreateAnnihilateOperatorMode[y];	
Sqrt[modeList[[operatorMode]]+1]*getPre[x]*getCreateAnnihilateOperatorPre[y]*ReplaceAt[Bra[getNum[x]],{modeList[[operatorMode]]->modeList[[operatorMode]]+1},{operatorMode}]
]

CenterDot[x_?isCreationOperator,y_?isDensity]:=
Module[
	{densityKet,operatorMode},
	densityKet=Ket[getNumDMLeft[y]];
	operatorMode=getCreateAnnihilateOperatorMode[x];
getPre[y]*getCreateAnnihilateOperatorPre[x]*SmallCircle[SuperDagger[Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]]\[CenterDot]densityKet,Bra[getNumDMRight[y]]]
]

CenterDot[x_?isAnnihilationOperator,y_?isDensity]:=
Module[
	{densityKet,operatorMode},
	densityKet=Ket[getNumDMLeft[y]];
	operatorMode=getCreateAnnihilateOperatorMode[x];
If[
	Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]\[CenterDot]densityKet=!=0,
	getPre[y]*getCreateAnnihilateOperatorPre[x]*SmallCircle[Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]\[CenterDot]densityKet,Bra[getNumDMRight[y]]],
	0
]
]

CenterDot[x_?isDensity,y_?isCreationOperator]:=
Module[
	{densityBra,operatorMode},
	densityBra=Bra[getNumDMRight[x]];
	operatorMode=getCreateAnnihilateOperatorMode[y];
If[
	densityBra\[CenterDot]SuperDagger[Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]]=!=0
	getPre[x]*getCreateAnnihilateOperatorPre[y]*SmallCircle[Ket[getNumDMLeft[x]],densityBra\[CenterDot]SuperDagger[Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]]],
	0
]
]

CenterDot[x_?isDensity,y_?isAnnihilationOperator]:=
Module[
	{densityBra,operatorMode},
	densityBra=Bra[getNumDMRight[x]];
	operatorMode=getCreateAnnihilateOperatorMode[y];
getPre[x]*getCreateAnnihilateOperatorPre[y]*SmallCircle[Ket[getNumDMLeft[x]],densityBra\[CenterDot]Subscript[\!\(\*OverscriptBox[\(a\), \(^\)]\), operatorMode]]
]



(*Sets power behaviour for creation and annihilation operators*)

(*CenterDot[x_?isCreationOperatorPower,y_?isKet]:=
CenterDot@@Join[ConstantArray[x/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[x]],{y}]

CenterDot[x_?isAnnihilationOperatorPower,y_?isKet]:=
CenterDot@@Join[ConstantArray[x/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[x]],{y}]

CenterDot[x_?isCreationOperatorPower,y_?isKet]:=
CenterDot@@Join[ConstantArray[x/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[x]],{y}]

CenterDot[x_?isAnnihilationOperatorPower,y_?isKet]:=
CenterDot@@Join[ConstantArray[x/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[x]],{y}]

CenterDot[x_?isBra,y_?isCreationOperatorPower]:=
CenterDot@@Join[{x},ConstantArray[y/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[y]]]

CenterDot[x_?isBra,y_?isAnnihilationOperatorPower]:=
CenterDot@@Join[{x},ConstantArray[y/.Power[a_,_]->a,getCreateAnnihilateOperatorPower[y]]]*)

(*TO DO: Power functionality for operators*)
(*TO DO: density matrices*)


(*Sets associativity.*)
CenterDot[x_?isOperatorComposition,y_?isKet]:=
Module[
	{args,len,tempProd},
	args=Join[List@@x,{y}];
	len=Length[args];
Do[
	tempProd=CenterDot[args[[-2]],args[[-1]]];
	args=Drop[args,-1];
	args[[-1]]=tempProd,
	{i,len-2}
];
Return[CenterDot[args[[-2]],args[[-1]]]]
]

CenterDot[x_?isBra,y__?isOperatorComposition]:=
Module[
	{args,len,tempProd},
	args=Join[{x},List@@y];
	len=Length[args];
Do[
	tempProd=CenterDot[args[[1]],args[[2]]];
	args=Drop[args,1];
	args[[1]]=tempProd,
	{i,len-2}
];
Return[CenterDot[args[[1]],args[[2]]]]
]

CenterDot[x__?isOperatorComposition,y_?isDensity]:=
Module[
	{args,len,tempProd},
	args=Join[List@@x,{y}];
	len=Length[args];
Do[
	tempProd=CenterDot[args[[-2]],args[[-1]]];
	args=Drop[args,-1];
	args[[-1]]=tempProd,
	{i,len-2}
];
Return[CenterDot[args[[-2]],args[[-1]]]]
]

CenterDot[x_?isDensity,y__?isOperatorComposition]:=
Module[
	{args,len,tempProd},
	args=Join[{x},List@@y];
	len=Length[args];
Do[
	tempProd=CenterDot[args[[1]],args[[2]]];
	args=Drop[args,1];
	args[[1]]=tempProd,
	{i,len-2}
];
Return[CenterDot[args[[1]],args[[2]]]]
]


(*Distributivity*)
CenterDot[x_?isOperator,y_?isOperatorSum]:=
Sum[
CenterDot[x,y[[i]]],
{i,Length[y]}
]

CenterDot[x_?isOperatorSum,y_?isOperator]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isOperatorSum,y_?isOperatorSum]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]

CenterDot[x_?isOperator,y_?isKetState]:=
Sum[
CenterDot[x,y[[i]]],
{i,Length[y]}
]

CenterDot[x_?isOperatorSum,y_?isKet]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isOperatorSum,y_?isKetState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]

CenterDot[x_?isBraState,y_?isOperator]:=
Sum[
CenterDot[x,y[[i]]],
{i,Length[y]}
]

CenterDot[x_?isBra,y_?isOperatorSum]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isBraState,y_?isOperatorSum]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]

CenterDot[x_?isOperator,y_?isDensityState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]

CenterDot[x_?isOperatorSum,y_?isDensity]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isOperatorSum,y_?isDensityState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]

CenterDot[x_?isDensity,y_?isOperatorSum]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]

CenterDot[x_?isDensityState,y_?isOperator]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isDensityState,y_?isOperatorSum]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]


(*Extend inner product to be associative for expectation-value type products, where the multiplication is of the form \[Psi]\[CenterDot]operator\[CenterDot]\[Rho].*)
CenterDot[x_?isBra|x_?isBraState|x_?isDensity|x_?isDensityState,z_?isOperator|z_?isOperatorSum,y_?isKet|y_?isKetState|y_?isDensity|y_?isDensityState]:=
CenterDot[x,CenterDot[z,y]]

(*Extend inner product to be associative for the special case of operator\[CenterDot]\[Rho]\[CenterDot]operator.*)
CenterDor[x_?isOperator|x_?isOperatorSum,z_?isDensity|z_?isDensityState,y_?isOperator|y_?isOperatorSum]:=
CenterDot[x,CenterDot[y,z]]


(*Extension to scalar behaviour for operators*)
CenterDot[x_?isOperator,Except[y_?isQuantum]]:=
y*x

CenterDot[Except[x_?isQuantum],y_?isOperator]:=
x*y

CenterDot[x_?isOperatorSum,Except[y_?isQuantum]]:=
Sum[
y*x[[i]],
{i,Length[x]}
]

CenterDot[Except[y_?isQuantum],y_?isOperatorSum]:=
Sum[
x*y[[j]],
{j,Length[y]}
]


(*Commutation relations*)


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


(*Define idempotency of the dagger operation.*)
SuperDagger[SuperDagger[x_]]:=
x


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


(*Trace out a density matrix according to the standard Fock eigenbasis. Quicker than generating a numerical matrix and using the inbuilt trace for large Hilbert spaces.*)
fullTrace[x_?isDensity]:=
If[
getNumDMLeft[x]===getNumDMRight[x],
getPre[x],
0
]

fullTrace[x_?isDensityState]:=
Sum[
fullTrace[x[[i]]],
{i,Length[x]}
]


(*xQuadratureMean[x_?isKet|x_?isDensity]:=
If[
	isKet[x],
	inputState=x\[CenterDot]x^\[Dagger]
];
fullTrace[

]

modeList=List[getNum[y]];
	operatorMode=getCreateAnnihilateOperatorMode[x];
Sqrt[modeList[[operatorMode]]+1]*getPre[y]*getCreateAnnihilateOperatorPre[x]*ReplaceAt[Ket[getNum[y]],{modeList[[operatorMode]]->modeList[[operatorMode]]+1},{operatorMode}]*)


Print["@MENTAT: All functions loaded. 
Welcome to MENTAT, a CAS system for Fock-state calculations in bosonic optical systems! Copyright N. Zaunders, University of Queensland, 2024. All rights reserved."]


End[];


EndPackage[];
