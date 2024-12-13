(* ::Package:: *)

(* ::Text:: *)
(*MENTAT package + documentation*)
(*Nov 2024*)
(*N. Zaunders, University of Queensland School of Mathematics and Physics*)
(*"MENTAT: Human Input, Machine Speeds"*)
(**)
(*----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------*)
(**)
(*To do:*)
(*- BUGFIX: why is a creation operator on the symbolic nth mode a valid form but an annihilation operator is not?*)
(*- Add logic for simplifying operator expressions like (3a . 3a . 3a) to 27(a.a.a) - not necessary, but makes things prettier. (Though this does occur sometimes??)*)


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
isCreationOperator::usage="Returns True if argument is a creation operator. Can have prefactors and be raised to any power."
isAnnihilationOperator::usage="Returns True if argument is an annihilation operator of the form. Can have prefactors and be raised to any power."
isCreationOperatorPower::usage="Returns True if argument is explicitly a creation operator raised to a power."
isAnnihilationOperatorPower::usage="Returns True if argument is explicitly an annihilation operator raised to a power."
isOperatorComposition::usage="Returns True is argument is a CenterDot composition of creation or annihilation operators."
isOperator::usage="Returns True if argument is a creation or annihilation operator or composition thereof."
isOperatorSum::usage="Returns True if argument is a sum of creation or annihilation operator raised to any power."
getCreateAnnihilateOperatorMode::usage="Returns mode index 'mode' of the input creation or annihilation operator \!\(\*SubscriptBox[OverscriptBox[\(a\), \(^\)], \(i\)]\)."
getCreateAnnihilateOperatorPower::usage="Returns power pow of the input creation or annihilation operator."
getCreateAnnihilateOperatorPre::usage="Returns prefactor pre of the input creation or annihilation operator."

CoherentState::usage="Takes complex amplitude \[Alpha] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the single-mode coherent state |\[Alpha]> up to cutoff. 
					  Optional argument allows cutoff to be chosen such that truncation error is below epsilon."
TwoModeSqueezedVacuumState::usage="Takes squeezing value \[Chi] and maximum Fock-state argument 'cutoff' and returns the ket state corresponding to the two-mode squeezed vacuum state |\[Chi]> up to cutoff.
							Optional argument allows cutoff to be chosen such that truncation error is below epsilon."
beamsplitter::usage="Implements beamsplitter functionality for a given ket-state x. Takes as input a single ket x, the beamsplitter transmissivity \[Tau], and the indices of the two modes to be mixed. Returns the ket state x' after the unitary beamsplitter is applied."
partialTrace::usage="Partial trace calculator. Takes as input a density matrix x of a state with n modes and a list tracedModesList of length < n describing the indices of the modes to be traced out. Returns a density matrix x' corresponding to x with the specified modes traced out."
generateNDArray::usage="Generates a list containing the basis Fock states of a Hilbert space corresponding to a system with number of modes numModes and maximum Fock-state argument cutoff. The list is of length (cutoff+1)^numModes."
matrixRepresent::usage="Takes a given density matrix x and the Hilbert space characteristics numModes, cutoff and produces a numeric matrix isomorphic to that density matrix in the basis space {|n1,n2,...,nNumModes>} for \!\(\*FormBox[\(\*SubscriptBox[\(n\), \(i\)] < cutoff\),
						TraditionalForm]\). Returns square matrix of dimension (cutoff+1)^numModes."
numericalMatrixFidelity::usage="Calculate the quantum fidelity of two positive semidefinite matrices x, y assuming they correspond to normalised and valid quantum states."
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
_Ket|Times[_,_Ket]
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
_Bra|Times[_,_Bra]
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
We use the SmallCircle operator \[SmallCircle] to fill the definition of a density state class.*)
isDensity[expr_]:=
MatchQ[
expr,
SmallCircle[_Ket,_Bra]|Times[_,SmallCircle[_Ket,_Bra]]
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


(*Returns the prefactor pre of expr when expr is a ket object,bra object,or density matrix object.
Different cases have to be put in to account for whether the object has a prefactor or not. We do this because using Mathematica's built-in Times to deal with scalar multiplication 
means objects multiplied by 1 or 0 will have a different code structure than those multiplied by an arbitrary quantity .*)
getPre[expr_?isKet|expr_?isBra|expr_?isDensity]:=
Which[
MatchQ[expr,_Ket],
	1,
MatchQ[expr,Times[_,_Ket]],
	expr/.Times[pre_,_Ket]->pre,
MatchQ[expr,_Bra],
	1,
MatchQ[expr,Times[_,_Bra]],
	expr/.Times[pre_,_Bra]->pre,
MatchQ[expr,SmallCircle[_Ket,_Bra]],
	1,
MatchQ[expr,Times[_,SmallCircle[_Ket,_Bra]]],
	expr/.Times[pre_,SmallCircle[_Ket,_Bra]]->pre
]


(*Retrieves the mode numbers num of the ket or bra state expr.  Returns num in List form.*)
getNum[expr_?isKet|expr_?isBra]:=
Which[
MatchQ[expr,_Ket],
	expr/.Ket[num_]->num,
MatchQ[expr,Times[_,_Ket]],
	expr/.Times[_,Ket[num_]]->num,
MatchQ[expr,__Bra],
	expr/.Bra[num_]->num,
MatchQ[expr,Times[_,_Bra]],
	expr/.Times[_,Bra[num_]]->num
]


(*Retrieves the mode numbers num of the ket state making up a density matrix object expr.  Returns num in List form.*)
getNumDMLeft[expr_?isDensity]:=
Which[
MatchQ[expr,SmallCircle[_Ket,_Bra]],
	expr/.SmallCircle[Ket[num_],_]->num,
MatchQ[expr,Times[_,SmallCircle[_Ket,_Bra]]],
	expr/.Times[_,SmallCircle[Ket[num_],_]]->num
]


(*Retrieves the mode numbers num of the bra state making up a density matrix object expr. Returns num in List form.*)
getNumDMRight[expr_?isDensity]:=
Which[
MatchQ[expr,SmallCircle[_Ket,_Bra]],
	expr/.SmallCircle[_,Bra[num_]]->num,
MatchQ[expr,Times[_,SmallCircle[_Ket,_Bra]]],
	expr/.Times[_,SmallCircle[_,Bra[num_]]]->num
]


(* ::Section:: *)
(*Inner product \[CenterDot]*)


(*Define operation of vector multiplication CenterDot for the vector inner product of two pure states. Returns scalar.
If x is a bra and y a ket, then return 0 if num of x != num of y, and return the product of the prefactors otherwise.*)
CenterDot[x_?isBra,y_?isKet]:=
If[
getNum[x]===getNum[y],
getPre[x]*getPre[y],
0
]


(*Define vector multiplication for vector outer product of two pure states. Returns density matrix.*)
CenterDot[x_?isKet,y_?isBra]:=
getPre[x]*getPre[y]*SmallCircle[Ket[getNum[x]],Bra[getNum[y]]]


(*Define vector multiplication for inner product between two density matrices. Returns density matrix.*)
CenterDot[x_?isDensity,y_?isDensity]:=
If[
getNumDMRight[x]===getNumDMLeft[y],
getPre[x]*getPre[y]*SmallCircle[Ket[getNumDMLeft[x]],Bra[getNumDMRight[y]]],
0
]


(*Define vector multiplication for inner product between left-bra and right-DM. Returns bra.*)
CenterDot[x_?isBra,y_?isDensity]:=
If[
getNum[x]===getNumDMLeft[y],
getPre[x]*getPre[y]*Bra[getNumDMRight[y]],
0
]


(*Define vector multiplication for inner product between left-DM and right-ket. Returns ket.*)
CenterDot[x_?isDensity,y_?isKet]:=
If[
getNumDMRight[x]===getNum[y],
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
Additional logic needs to be put in for the degenerate case where a single ket or bra is multplied by a ket state or bra state. 
This probably isn't very efficient, since the amount of different cases grows exponentially. Fortunately we only really need to consider bras, kets, and density matrices. 
Another important reason for doing it this way is that we can explicitly exclude unwanted scenarios, like ket\[CenterDot]DM or DM\[CenterDot]bra.*)
CenterDot[x_?isBraState,y_?isKetState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isBra,y_?isKetState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isBraState,y_?isKet]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isKetState,y_?isBraState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isKet,y_?isBraState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isKetState,y_?isBra]:=
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
CenterDot[x_?isBra,y_?isDensityState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isBraState,y_?isDensity]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isDensityState,y_?isKetState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isDensity,y_?isKetState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isDensityState,y_?isKet]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isDensityState,y_?isDensityState]:=
Sum[
CenterDot[x[[i]],y[[j]]],
{i,Length[x]},
{j,Length[y]}
]
CenterDot[x_?isDensity,y_?isDensityState]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
]
CenterDot[x_?isDensityState,y_?isDensity]:=
Sum[
CenterDot[x[[i]],y],
{i,Length[x]}
]


(*Extend inner product to be associative for expectation-value type products, where the multiplication is of the form \[LeftAngleBracket]a|\[CenterDot]b\[CenterDot]|c\[RightAngleBracket].*)
CenterDot[x_?isBra|x_?isBraState,z_,y_?isKet|y_?isKetState]:=
CenterDot[x,CenterDot[z,y]]


(* ::Section:: *)
(*Outer product \[CircleTimes]*)


(*Base definition of the tensor product. *)
CircleTimes[x_?isKet,y_?isKet]:=
getPre[x]*getPre[y]*Ket[Join[getNum[x],getNum[y]]]

CircleTimes[x_?isBra,y_?isBra]:=
getPre[x]*getPre[y]*Bra[Join[getNum[x],getNum[y]]]

CircleTimes[x_?isDensity,y_?isDensity]:=
getPre[x]*getPre[y]*SmallCircle[Ket[Join[getNumDMLeft[x],getNumDMLeft[y]]],Bra[Join[getNumDMRight[x],getNumDMRight[y]]]]


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
	{args,len},
	args=List[x];
	len=Length[args];
Do[
	args[[-2]]=CircleTimes[args[[-2]],args[[-1]]];
	args=Drop[args,-1],
	{i,len-2}
];
CircleTimes[args[[-2]],args[[-1]]]
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


(*Define identifier function for creation operators. Returns True if expr is a creation operator.
NOTE: We use HoldPattern on the pattern to prevent recursion happening via intersection with our definition of SuperDagger - if we want to define SuperDagger behaviour on e.g.
creation operators, then using isCreationOperator in the pattern specification calls SuperDagger, which calls isCreationOperator, which calls SuperDagger, etc. Unevaluated makes 
the pattern matching literal, i.e. MatchQ only returns true if the input form is explicitly and literally the one given.*)
isCreationOperator[expr_]:=
MatchQ[
	expr,
	HoldPattern[
		 SuperDagger[Subscript[OverHat[a_Symbol],_Integer]]
		|Times[_,SuperDagger[Subscript[OverHat[a_Symbol],_Integer]]]
	]
]

(*Define identifier function for powers of annihilation operators. Returns True if expr is an annihilation operator.*)
isAnnihilationOperator[expr_]:=
MatchQ[
	expr,
	HoldPattern[
		 Subscript[OverHat[a_Symbol],_Integer]
		|Times[_,Subscript[OverHat[a_Symbol],_Integer]]
	]
]

(*Define identifier function for powers of annihilation operators. Returns True if expr is an annihilation operator.*)
isCreationOperatorPower[expr_]:=
MatchQ[
	expr,
	HoldPattern[
		Power[SuperDagger[Subscript[OverHat[_],_]],_]
		|Times[_,Power[SuperDagger[Subscript[OverHat[_],_]],_]]
	]
]

isAnnihilationOperatorPower[expr_]:=
MatchQ[
	expr,
	HoldPattern[
		Power[Subscript[OverHat[_],_],_]
		|Times[_,Power[Subscript[OverHat[_],_],_]]
	]
]

(*Define identifier function for a composition of Fock state operators, i.e. Subscript[Overscript[a, ^], i]^\[Dagger]\[CenterDot]Subscript[Overscript[a, ^], j]^\[Dagger]\[CenterDot]..., etc. Returns True if expr is an
arbitrary series of creation or annihilation operators CenterDot'd with each other.
Works by checking truth of isCreation/isAnnihilation for each element in the Sequence passed to CenterDot. The Sequence of booleans is then checked via And; if every element returns true, 
then it must be a composition of operators.*)
isOperatorComposition[expr_]:=
MatchQ[
	expr,
	CenterDot[x__/;(And@@(MatchQ[_?isCreationOperator|_?isAnnihilationOperator|_?isCreationOperatorPower|_?isAnnihilationOperatorPower]/@(List[x])))]
]

(*Define identifier function for any type of Fock state operator. Returns True if expr is a creation or annihilation operator or power or composition thereof.*)
isOperator[expr_]:=
MatchQ[
	expr,
	_?isCreationOperator|_?isAnnihilationOperator|_?isCreationOperatorPower|_?isAnnihilationOperatorPower|_?isOperatorComposition
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

(*Returns the mode index that the creation or annihilation operator x should operate on.
Note that HoldPattern is required here, and in other functions of similar type, wherever pattern matching is performed on a pattern using SuperDagger to prevent recursion issues.*)
getCreateAnnihilateOperatorMode[expr_?isCreationOperator|expr_?isCreationOperatorPower|expr_?isAnnihilationOperator|expr_?isAnnihilationOperatorPower]:=
Which[
isCreationOperator[expr],
	Which[
		MatchQ[expr,HoldPattern[SuperDagger[Subscript[OverHat[_],mode_Integer]]]],
			expr/.HoldPattern[SuperDagger[Subscript[OverHat[_],mode_Integer]]]->mode,
		MatchQ[expr,HoldPattern[Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]]],
			expr/.HoldPattern[Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]]->mode
	],
isCreationOperatorPower[expr],
	Which[		
		MatchQ[expr,HoldPattern[Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]],
			expr/.HoldPattern[Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]->mode,
		MatchQ[expr,HoldPattern[Times[_,Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]]],
			expr/.HoldPattern[Times[_,Power[SuperDagger[Subscript[OverHat[_],mode_Integer]],_]]]->mode
	],
isAnnihilationOperator[expr],
	Which[
		MatchQ[expr,Subscript[OverHat[_],mode_Integer]],
			expr/.Subscript[OverHat[_],mode_Integer]->mode,	
		MatchQ[expr,Times[_,Subscript[OverHat[_],mode_Integer]]],
			expr/.Times[_,Subscript[OverHat[_],mode_Integer]]->mode
	],
isAnnihilationOperatorPower[expr],
	Which[
		MatchQ[expr,Power[Subscript[OverHat[_],mode_Integer],_]],
			expr/.Power[Subscript[OverHat[_],mode_Integer],_]->mode,	
		MatchQ[expr,Times[_,Power[Subscript[OverHat[_],mode_Integer],_]]],
			expr/.Times[_,Power[Subscript[OverHat[_],mode_Integer],_]]->mode
	]
]

(*Returns the power that the creation or annihilation operator expr is raised to.*)
getCreateAnnihilateOperatorPower[expr_?isCreationOperator|expr_?isCreationOperatorPower|expr_?isAnnihilationOperator|expr_?isAnnihilationOperatorPower]:=
Which[
isCreationOperator[expr],
	Which[
	MatchQ[expr,HoldPattern[SuperDagger[Subscript[OverHat[_],mode_Integer]]]],
		1,
	MatchQ[expr,HoldPattern[Times[_,SuperDagger[Subscript[OverHat[_],mode_Integer]]]]],
		1
	],
isCreationOperatorPower[expr],
	Which[
	MatchQ[expr,HoldPattern[Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]],
		expr/.HoldPattern[Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]->pow,
	MatchQ[expr,HoldPattern[Times[_,Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]]],
		expr/.HoldPattern[Times[_,Power[SuperDagger[Subscript[OverHat[_],_]],pow_Integer]]]->pow
	],
isAnnihilationOperator[expr],
	Which[
	MatchQ[expr,Subscript[OverHat[_],_]],
		1,
	MatchQ[expr,Times[_,Subscript[OverHat[_],_]]],
		1
	],
isAnnihilationOperatorPower[expr],
	Which[
	MatchQ[expr,Power[Subscript[OverHat[_],_],pow_Integer]],
		expr/.Power[Subscript[OverHat[_],_],pow_Integer]->pow,
	MatchQ[expr,Times[_,Power[Subscript[OverHat[_],_],pow_Integer]]],
		expr/.Times[_,Power[Subscript[OverHat[_],_],pow_Integer]]->pow
	]
]

(*Returns the prefactor of the creation or annihilation operator.*)
getCreateAnnihilateOperatorPre[expr_?isCreationOperator|expr_?isCreationOperatorPower|expr_?isAnnihilationOperator|expr_?isAnnihilationOperatorPower]:=
Which[
isCreationOperator[expr],
	Which[
	MatchQ[expr,HoldPattern[SuperDagger[Subscript[OverHat[_],_]]]],
		1,
	MatchQ[expr,HoldPattern[Times[pre_,SuperDagger[Subscript[OverHat[_],_]]]]],
		expr/.HoldPattern[Times[pre_,SuperDagger[Subscript[OverHat[_],_]]]]->pre
	],
isCreationOperatorPower[expr],
	Which[
	MatchQ[expr,HoldPattern[Power[SuperDagger[Subscript[OverHat[_],_]],_]]],
		1,
	MatchQ[expr,HoldPattern[Times[pre_,Power[SuperDagger[Subscript[OverHat[_],_]],_]]]],
		expr/.HoldPattern[Times[pre_,Power[SuperDagger[Subscript[OverHat[_],_]],_]]]->pre
	],
isAnnihilationOperator[expr],
	Which[
	MatchQ[expr,Subscript[OverHat[_],_]],
		1,
	MatchQ[expr,Times[pre_,Subscript[OverHat[_],_]]],
		expr/.Times[pre_,Subscript[OverHat[_],_]]->pre
	],
isAnnihilationOperatorPower[expr],
	Which[
	MatchQ[expr,Power[Subscript[OverHat[_],_],_]],
		1,
	MatchQ[expr,Times[pre_,Power[Subscript[OverHat[_],_],_]]],
		expr/.Times[pre_,Power[Subscript[OverHat[_],_],_]]->pre
	]
]


(*Define raising and lowering operations. We re-use the operator \[CenterDot] to signify interaction with an operator.
To give multimodal functionality, each operator is equipped with an integer subscript which corresponds to the mode it should operate on.*)
CenterDot[x_?isCreationOperator,y_?isKet]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[x], 
	modeList=getNum[y]
},
modeList[[operatorMode]]=modeList[[operatorMode]]+1;
Sqrt[modeList[[operatorMode]]]*getPre[y]*getCreateAnnihilateOperatorPre[x]*Ket[modeList]
]

CenterDot[x_?isAnnihilationOperator,y_?isKet]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[x],
	modeList=getNum[y]
},
modeList[[operatorMode]]=modeList[[operatorMode]]-1;
If[
	modeList[[operatorMode]]=!=-1,
	Sqrt[modeList[[operatorMode]]+1]*getPre[y]*getCreateAnnihilateOperatorPre[x]*Ket[modeList],
	0
]
]

CenterDot[x_?isBra,y_?isCreationOperator]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[y],
	modeList=getNum[x]
},
modeList[[operatorMode]]=modeList[[operatorMode]]-1;
If[
	modeList[[operatorMode]]=!=-1,
	Sqrt[modeList[[operatorMode]]+1]*getPre[x]*getCreateAnnihilateOperatorPre[y]*Bra[modeList],
	0
]
]

CenterDot[x_?isBra,y_?isAnnihilationOperator]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[y], 
	modeList=getNum[x]
},
modeList[[operatorMode]]=modeList[[operatorMode]]+1;
Sqrt[modeList[[operatorMode]]]*getPre[x]*getCreateAnnihilateOperatorPre[y]*Bra[modeList]
]

CenterDot[x_?isCreationOperator,y_?isDensity]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[x],
	ketModeList=getNumDMLeft[y]
},
ketModeList[[operatorMode]]=ketModeList[[operatorMode]]+1;
Sqrt[ketModeList[[operatorMode]]]*getPre[y]*getCreateAnnihilateOperatorPre[x]*SmallCircle[Ket[ketModeList],Bra[getNumDMRight[y]]
]
]

CenterDot[x_?isAnnihilationOperator,y_?isDensity]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[x],
	ketModeList=getNumDMLeft[y]
},
ketModeList[[operatorMode]]=ketModeList[[operatorMode]]-1;
If[
	ketModeList[[operatorMode]]=!=-1,
	Sqrt[ketModeList[[operatorMode]]+1]*getPre[y]*getCreateAnnihilateOperatorPre[x]*SmallCircle[Ket[ketModeList],Bra[getNumDMRight[y]]],
	0
]
]

CenterDot[x_?isDensity,y_?isCreationOperator]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[y],
	braModeList=getNumDMRight[x]
},
braModeList[[operatorMode]]=braModeList[[operatorMode]]-1;
If[
	braModeList[[operatorMode]]=!=-1,
	Sqrt[braModeList[[operatorMode]]+1]*getPre[x]*getCreateAnnihilateOperatorPre[y]*SmallCircle[Ket[getNumDMLeft[x]],Bra[braModeList]],
	0
]
]

CenterDot[x_?isDensity,y_?isAnnihilationOperator]:=
Module[{
	operatorMode=getCreateAnnihilateOperatorMode[y],
	braModeList=getNumDMRight[x]
},
braModeList[[operatorMode]]=braModeList[[operatorMode]]+1;
Sqrt[braModeList[[operatorMode]]]*getPre[x]*getCreateAnnihilateOperatorPre[y]*SmallCircle[Ket[getNumDMLeft[x]],Bra[braModeList]]
]



(*Sets power behaviour for creation and annihilation operators. *)
CenterDot[x_?isCreationOperatorPower,y_?isKet]:=
getCreateAnnihilateOperatorPre[x]*CenterDot@@Join[
	ConstantArray[SuperDagger[Subscript[OverHat[a],getCreateAnnihilateOperatorMode[x]]],getCreateAnnihilateOperatorPower[x]],
	{y}
]

CenterDot[x_?isAnnihilationOperatorPower,y_?isKet]:=
getCreateAnnihilateOperatorPre[x]*CenterDot@@Join[
	ConstantArray[Subscript[OverHat[a],getCreateAnnihilateOperatorMode[x]],getCreateAnnihilateOperatorPower[x]],
	{y}
]

CenterDot[x_?isBra,y_?isCreationOperatorPower]:=
getCreateAnnihilateOperatorPre[y]*CenterDot@@Join[
	{x},
	ConstantArray[SuperDagger[Subscript[OverHat[a],getCreateAnnihilateOperatorMode[y]]],getCreateAnnihilateOperatorPower[y]]
]

CenterDot[x_?isBra,y_?isAnnihilationOperatorPower]:=
getCreateAnnihilateOperatorPre[y]*CenterDot@@Join[
	{x},
	ConstantArray[Subscript[OverHat[a],getCreateAnnihilateOperatorMode[y]],getCreateAnnihilateOperatorPower[y]]
]


(*Sets associativity. For a given composition of N operators followed by a Ket object, the last operator in the sequence is applied to the Ket and the second-to-last item in the total
sequence is updated with the result. The last item in the total sequence (the original Ket) is then deleted, leaving a composition of N-1 operators followed by a Ket. This repeats
until there is only a single operator left, at which point the function returns the product of the operator applied to the Ket.*)
CenterDot[x_?isOperatorComposition,y_?isKet|y_?isDensity]:=
Module[{
	args=Join[List@@x,{y}],
	len
},
len=Length[args];
Do[
	args[[-2]]=CenterDot[args[[-2]],args[[-1]]];
	args=Drop[args,-1],
	{i,len-2}
];
Return[CenterDot[args[[-2]],args[[-1]]]]
]

(*This special case is needed for when the operator part is not enclosed in brackets.*)
CenterDot[x__/;isOperatorComposition[CenterDot[x]],y_?isKet|y_?isDensity]:=
CenterDot[CenterDot[x],y]


CenterDot[x_?isBra|x_?isDensity,y_?isOperatorComposition]:=
Module[{
	args=Join[{x},List@@y],
	len
},
len=Length[args];
Do[
	args[[2]]=CenterDot[args[[1]],args[[2]]];
	args=Drop[args,1],
	{i,len-2}
];
Return[CenterDot[args[[1]],args[[2]]]]
]

CenterDot[x_?isBra|x_?isDensity,y__/;isOperatorComposition[CenterDot[y]]]:=
CenterDot[x,CenterDot[y]]


(*Simplify compositions of compositions and compositions of sums.*)
CenterDot[x_?isOperatorComposition,y_?isOperatorComposition]:=
CenterDot@@Join[]


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
CenterDot[x[[i]],y],
{i,Length[x]}
]

CenterDot[x_?isBra,y_?isOperatorSum]:=
Sum[
CenterDot[x,y[[j]]],
{j,Length[y]}
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
CenterDot[x_?isOperator|x_?isOperatorSum,z_?isDensity|z_?isDensityState,y_?isOperator|y_?isOperatorSum]:=
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

CenterDot[Except[x_?isQuantum],y_?isOperatorSum]:=
Sum[
x*y[[j]],
{j,Length[y]}
]


(* ::Subsection:: *)
(*Conjugate transpose (dagger)*)


(*Define idempotency of the dagger operation.*)
SuperDagger[SuperDagger[x_]]:=
x


(*Defines operation of the conjugate tranpose operation ^Dagger[*]. For an input quantum state x, conjugate the prefactor and transpose the quantum state appropriately.*) 
SuperDagger[x_?isKet]:=
Times[Conjugate[getPre[x]],Bra[getNum[x]]]

SuperDagger[x_?isBra]:=
Conjugate[getPre[x]]*Ket[getNum[x]]

SuperDagger[x_?isDensity]:=
Conjugate[getPre[x]]*SmallCircle[Ket[getNumDMRight[x]],Bra[getNumDMLeft[x]]]


(*Define distributivity of the conjugate transpose for quantum states.*)
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


(*Defines operation of the conjugate tranpose for creation operators.*) 
SuperDagger[SuperDagger[Subscript[OverHat[a_Symbol],mode_Integer]]]:=
Subscript[OverHat[a],mode]

SuperDagger[Times[pre_,SuperDagger[Subscript[OverHat[a_Symbol],mode_Integer]]]]:=
Times[Conjugate[pre],Subscript[OverHat[a],mode]]

SuperDagger[Power[SuperDagger[Subscript[OverHat[a_Symbol],mode_Integer]],pow_Integer]]:=
Power[Subscript[OverHat[a],mode],pow]

SuperDagger[Times[pre_,Power[SuperDagger[Subscript[OverHat[a_Symbol],mode_Integer]],pow_Integer]]]:=
Times[Conjugate[pre],Power[Subscript[OverHat[a],mode],pow]]


(*Defines operation of the conjugate tranpose for annihilation operators.
By the way operators are defined in MENTAT, we don't need to define a behaviour for the base case - Mathematica handles it automatically.*) 
SuperDagger[Times[pre_,Subscript[OverHat[a_Symbol],mode_Integer]]]:=
Times[Conjugate[pre],SuperDagger[Subscript[OverHat[a],mode]]]

SuperDagger[Power[Subscript[OverHat[a_Symbol],mode_Integer],pow_Integer]]:=
Power[SuperDagger[Subscript[OverHat[a],mode]],pow]

SuperDagger[Times[pre_,Power[Subscript[OverHat[a_Symbol],mode_Integer],pow_Integer]]]:=
Times[Conjugate[pre],Power[SuperDagger[Subscript[OverHat[a],mode]],pow]]


(*Define distributivity of the conjugate transpose for operators over addition.*)
SuperDagger[x_?isOperatorSum]:=
Sum[
	SuperDagger[x[[i]]],
	{i,Length[x]}
]


(*Define distributivity of the conjugate transpose for operators over composition.*)
SuperDagger[x_?isOperatorComposition]:=
CenterDot@@(SuperDagger/@(List@@x))


(*Conjugate transpose behaviour for scalars, i.e. reduction to complex conjugation.*)
SuperDagger[Except[x_?isQuantum]]:=
Conjugate[x]


(* ::Section:: *)
(*Functionality*)


(*Returns a ket state describing the coherent state Ket[\[Alpha]]for a complex amplitude \[Alpha].
Optional arguments allow the returned state to be specified either in terms of an explicit maximum Fock-state cutoff n, or alternatively a fraction F < 1 may be specified and
the returned state will be the smallest-dimension state with norm >= F.*)
Options[CoherentState]={n->-1,F->-1};

CoherentState[\[Alpha]_,OptionsPattern[]]:=
Module[{
	norm=0.0,
	i=0,
	state
},
Which[
	OptionValue[n]!=-1,
		Sum[Exp[-(Abs[\[Alpha]]^2/2)] \[Alpha]^k/\[Sqrt](k!) Ket[{k}],{k,0,OptionValue[n]}],
		
	OptionValue[F]!=-1&&NumericQ[\[Alpha]],
		While[norm<=OptionValue[F],
			state=Sum[Exp[-(Abs[\[Alpha]]^2/2)] \[Alpha]^k/\[Sqrt](k!) Ket[{k}],{k,0,i}];
			norm=SuperDagger[state]\[CenterDot]state;
			i=i+1;
		];
		Return[state]
]
]


(*Returns a two-mode squeezed vacuum state Ket[\[Chi]] for squeezing \[Chi] = tanh(r) and maximum Fock-state element cutoff.*)
Options[TwoModeSqueezedVacuumState]={n->-1,F->-1};

TwoModeSqueezedVacuumState[\[Chi]_,OptionsPattern[]]:=
Module[{
	norm=0.0,
	i=0,
	state
},
Which[
	OptionValue[n]!=-1,
		Sum[\[Sqrt](1-\[Chi]^2) \[Chi]^k Ket[{k,k}],{k,0,OptionValue[n]}],
		
	OptionValue[F]!=-1&&NumericQ[\[Chi]],
		While[norm<=OptionValue[F],
			state=Sum[\[Sqrt](1-\[Chi]^2) \[Chi]^k Ket[{k,k}],{k,0,i}];
			norm=SuperDagger[state]\[CenterDot]state;
			i=i+1;
		];
		Return[state]
]
]


(*Implements beamsplitter functionality for a given ket-state x. Takes as input a single ket x, the beamsplitter transmissivity \[Tau], and a list mixList containing the indices of the two modes
to be mixed. Returns the ket state x' after the unitary beamsplitter is applied.*)
beamsplitter[x_?isKet,mixList_List,\[Tau]_]:=
Module[{
	modeList=getNum[x],
	n=getNum[x][[mixList[[1]]]],
	m=getNum[x][[mixList[[2]]]]
},
Expand[
	Sum[
		modeList[[mixList[[1]]]]=n-i+m-j;
		modeList[[mixList[[2]]]]=i+j;
		getPre[x]*(-1)^i*(Sqrt[(i+j)!]*Sqrt[m!]*Sqrt[n!]*Sqrt[(n-i+m-j)!])/(i!*j!*(m-j)!*(n-i)!)*Sqrt[(\[Tau])^(n-i+j)*(1-\[Tau])^(i+m-j)]*Ket[modeList],
		{i,0,n},
		{j,0,m}
	]
]
]


(*Defines distributivity for the beamsplitter operation.*)
beamsplitter[x_?isKetState,mixList_List,\[Tau]_]:=
Sum[
	beamsplitter[x[[i]],mixList,\[Tau]],
	{i,Length[x]}
]


(*Partial trace calculator. Takes as input a density matrix x of a state with n modes and a list tracedModesList of length < n describing the indices of the modes to be traced out. Returns a density matrix x' corresponding to x with the specified modes traced out.
Note: this is less of a mathematically rigorous definition and more of a computational one. Faster and easier, but it's a shame to lose the mathematical backing.*)
partialTrace[x_?isDensity,tracedModesList_List]:=
Module[{
	numLeft=getNumDMLeft[x],
	numRight=getNumDMRight[x],
	tracedNumLeft,
	tracedNumRight
},
tracedNumLeft=Delete[numLeft,List/@tracedModesList];
tracedNumRight=Delete[numRight,List/@tracedModesList];
If[
	numLeft[[tracedModesList]]===numRight[[tracedModesList]],
	getPre[x]*SmallCircle[Ket[tracedNumLeft],Bra[tracedNumRight]],
	0
]
]


(*Define distributivity of the partial trace.*)
partialTrace[x_?isDensityState,tracedModesList_List]:=
Sum[
	partialTrace[x[[i]],tracedModesList],
	{i,Length[x]}
]


(*Generates a list containing the basis Fock states of a Hilbert space corresponding to a system with number of modes numModes and maximum Fock-state argument cutoff.
The list is of length (cutoff+1)^numModes.*)
generateNDArray[numModes_,cutoff_]:=
Module[{
	basisVectorSum=Sum[Ket[{i}],{i,0,cutoff}]
},
Do[
	basisVectorSum=(basisVectorSum)\[CircleTimes](Sum[Ket[{i}],{i,0,cutoff}]),
	{j,1,numModes-1}
];
Return[List@@basisVectorSum]
]


(*Takes a given density matrix x and the Hilbert space characteristics numModes, cutoff and produces a numeric matrix isomorphic to that density matrix in the basis space {|n1,n2,...,nNumModes>} for Subscript[n, i]<cutoff. Returns square matrix of dimension (cutoff+1)^numModes.*)
matrixRepresent[x_?isDensityState,numModes_,cutoff_]:=
Module[{
	rhoMatrix=ConstantArray[0,{(cutoff+1)^numModes,(cutoff+1)^numModes}],
	basisList=generateNDArray[numModes,cutoff]
},
Do[
	rhoMatrix[[i,j]]=SuperDagger[basisList[[i]]]\[CenterDot]x\[CenterDot]basisList[[j]],
	{i,1,(cutoff+1)^numModes},
	{j,1,(cutoff+1)^numModes}
];
Return[rhoMatrix]
]


(*Calculate the quantum fidelity of two positive semidefinite matrices x, y assuming they correspond to normalised and valid quantum states.*)
numericalMatrixFidelity[x_,y_]:=
Tr[MatrixPower[MatrixPower[x,1/2] . y . MatrixPower[x,1/2],1/2]]^2


(*Calculate the von Neumann entropy of a given density matrix state.*)
vonNeumannEntropy[x_?isDensityState,numModes_,cutoff_]:=
Module[{
	numericalMatrix=matrixRepresent[x,numModes,cutoff],
	eigenvalues
},
eigenvalues=Eigenvalues[numericalMatrix];
Sum[
	If[
		eigenvalues[[i]]===0,
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


Print["@MENTAT: All functions loaded. 
Welcome to MENTAT, a computer algebra system for Fock-state calculations in bosonic optical systems! Copyright N. Zaunders, University of Queensland, 2024. All rights reserved."]


End[];


EndPackage[];
