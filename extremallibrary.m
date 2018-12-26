BeginPackage["extremallibrary`"]
Get["./Quantum.m"]

Coherent::usage = "Defines a finite approximation to a coherent state. Requires a complex number as first input and an
integer for the dimension to use."
PhaseNumberOperatorMatrix::usage = "Induces a phase depending in the number of the Fock space. Requires the phase angle and
the dimension of the state it is acting on."
PsiTheta::usage = "Gives the number operator exponential applied to the coherent state. Requires the same as Coherent plus the phase."
PsiThetaNormalized::usage = "Normalized version of the previous state."
PsiThetaDM::usage = "Density matrix of the normalized state."
expectationN::usage = "Expectation value of the number of photons of the field in terms of the frequency and temperature."
Thermal::usage = "Thermal state that depends on the dimension that one is approximating to, the field frequency and the temperature."
ThermalMod::usage = "Modifiel Thermal state by the phase operation."
kBoltzmann::usage = "Boltzmann constant."
hBar::usage = "Value of the Planck constant divided by 2Pi."
QubitState::usage = "The state of a Qubit given by two angles in the Bloch sphere: Eta and Theta."
QBDM::usage = "Density matrix correspondent to the Qubit."
CoherentPlusThermalState::usage = "The convex sum of a Coherent plus a thermal state modified by the phase operator."
GellMann::usage = "Yields a table with all the Gell-Mann Matrices of a given dimension. You must input the dimension. This
function was found in the blog: https://blog.suchideas.com/2016/04/sun-gell-mann-matrices-in-
mathematica/ ."
CanonicalProjector::usage = "Canonical basis for the matrix space of dimension dim."
POVM::usage = "Generates a random POVM of dimension HilbertDim with Outcomedim the number of outcomes."
Gammapdf::usage = "We use the Gamma distribution for our calculation."
Sobra::usage = "The Fisher Information correspondent to the a priori distribution."
eivect::usage = "Eigenvector j of unPOVM element i."
eivalues::usage = "Eigenvalue j of unPOVM element i."
projector::usage = "Projector to eigenvector i,j."
POVMUnit::usage = "The normalized version of POVM."
PderivCoherentPlusTher::usage = "Derivative of the probability distribution given the extremal POVM for Coherent plus Thermal case. Requires the Phase, the number of the POVM element and a complex number for the Coherent state."
FisherCoherentPlusTher::usage = "Fisher Information given the extremal POVM for the Coherent plus Thermal case. Requires the Phase and the initial Complex number of the coherent state."
PderivQubit::usage = "Derivative of the probability distribution given the extremal POVM for Qubit case. 
Requires two angles in the Bloch sphere."
FisherQubit::usage = "Fisher Information given the extremal POVM for the Qubit case. Requires two angles."
Gaussianpdf::usage = "The Gaussian distribution probability density function."
NormCnsnt::usage = "Normalization constant for Gaussianpdf."
SobraGaussian::usage = "Fisher information for a priori Gaussian distribution."

Begin["Private`"]      
(*Done {{{*)
kBoltzmann = 1.38*^-23;
hBar = (6.626*^-34)/2 \[Pi];

Coherent[ComplexCoherent_, Dimension_] := Exp[-(1/2)*((Norm[ComplexCoherent])^2)] Sum[((ComplexCoherent^(n - 1))/Sqrt[((n - 1)!)])*UnitVector[Dimension + 1, n], {n, 1, Dimension + 1}];
PhaseNumberOperatorMatrix[Theta_, Dimension_] := Table[Chop[KroneckerDelta[i, j]*Exp[I*Theta*(i - 1)]], {i, Dimension + 1}, {j, Dimension + 1}];
PsiTheta[ComplexCoherent_, Theta_, Dimension_] := PhaseNumberOperatorMatrix[Theta, Dimension].Coherent[ComplexCoherent, Dimension];
PsiThetaNormalized[ComplexCoherent_, Theta_, Dimension_] := PsiTheta[ComplexCoherent, Theta, Dimension]/(PsiTheta[Conjugate[ComplexCoherent], -Theta, Dimension].PsiTheta[ComplexCoherent, Theta, Dimension]);
PsiThetaDM[ComplexCoherent_, Theta_, Dimension_] := TensorProduct[PsiThetaNormalized[ComplexCoherent, Theta, Dimension],PsiThetaNormalized[Conjugate[ComplexCoherent], -Theta, Dimension]];
expectationN[FieldFrequency_,  Temperature_] := 1/(Exp[(hBar*FieldFrequency)/(kBoltzmann*Temperature)] - 1);
Thermal[Dimension_, FieldFrequency_, Temperature_] := Sum[TensorProduct[UnitVector[Dimension + 1, n],UnitVector[Dimension + 1, n]]*expectationN[FieldFrequency, Temperature]^n/(1 + expectationN[FieldFrequency, Temperature])^(n + 1), {n, Dimension}];
ThermalMod[Dimension_, FieldFrequency_, Temperature_, Theta_] := (PhaseNumberOperatorMatrix[-Theta, Dimension].Thermal[Dimension, FieldFrequency, Temperature]).PhaseNumberOperatorMatrix[Theta, Dimension];
QubitState[Eta_, Theta_] := Exp[-I*Theta/2]*Cos[Eta/2]*UnitVector[2, 1] + Exp[I*Theta/2]*Sin[Eta/2]*UnitVector[2, 2];
QBDM[Eta_, Theta_] := TensorProduct[QubitState[Eta, Theta],QubitState[Eta, -Theta]];
CoherentPlusThermalState[Dimension_, Theta_, ComplexCoherent_, FieldFrequency_, Temperature_, MixConstant_] := MixConstant*PsiThetaDM[ComplexCoherent, Theta, Dimension] + (1 - MixConstant)*ThermalMod[Dimension, FieldFrequency, Temperature, Theta];
GellMann[Dimension_] := Flatten[Table[(*Symmetric case*)
    SparseArray[{{j, k} -> 1, {k, j} -> 1}, {Dimension, Dimension}], {k, 2, Dimension}, {j,  1, k - 1}], 1]~Join~Flatten[
   Table[(*Antisymmetric case*)
    SparseArray[{{j, k} -> -I, {k, j} -> +I}, {Dimension, Dimension}], {k, 2, Dimension}, {j, 1, k - 1}], 1]~Join~
    Table[(*Diagonal case*)
      Sqrt[2/l/(l + 1)] SparseArray[Table[{j, j} -> 1, {j, 1, l}]~Join~{{l + 1, l + 1} -> -l}, {Dimension, Dimension}], {l, 1, Dimension - 1}];
CanonicalProjector[iterator_,Outcomedim_] := TensorProduct[UnitVector[Outcomedim, iterator],UnitVector[Outcomedim, iterator]];
POVM[Iterator_,HilbertDim_,Outcomedim_] := If[IntegerQ[Iterator] == True,  
    HRandomMatrix = CUEMember[Outcomedim*(HilbertDim)]; 
    X = {};
      Y = {};
      For[s = 0, s < HilbertDim, s++,
        Y = {};
        For[j = 0, j < HilbertDim, j++,
          AppendTo[Y, 
            Chop[Sum[(Transpose[
            Conjugate[HRandomMatrix]].(KroneckerProduct[IdentityMatrix[HilbertDim], 
                          CanonicalProjector[Iterator,Outcomedim]]).HRandomMatrix)[[i + s*Outcomedim, i + j*Outcomedim]], {i, 1, Outcomedim}]]]];
        AppendTo[X, Y];]; Return[X/Outcomedim]];
Gammapdf[Alpha_, Beta_, x_] := Exp[-x/Beta]*x^(Alpha - 1)*Beta^(-Alpha)/Gamma[Alpha];
Sobra = NIntegrate[Power[D[Gammapdf[4, 1.5, x], x], 2]/Gammapdf[4, 1.5, x], {x, 0, Infinity}];
eivect[unPOVM_,i_, j_] := Eigenvectors[unPOVM[[i]]][[j]];
eivalues[unPOVM_,i_, j_] := Eigenvalues[unPOVM[[i]]][[j]];
projector[unPOVM_,i_, j_] := TensorProduct[eivect[unPOVM,i, j],Conjugate[eivect[unPOVM,i, j]]];
POVMUnit[unPOVM_,Iterator_] := unPOVM[[Iterator]]/Tr[unPOVM[[Iterator]]];
PderivCoherentPlusTher[Extremal_,ThetaPhase_, Iterator_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_,MixConstant_] := 
Module[{DerivativeVar},D[Tr[CoherentPlusThermalState[HilbertDim -1, DerivativeVar, ComplexCoherent, FieldFrequency, Temperature, MixConstant].Extremal[[Iterator]]], DerivativeVar] /.DerivativeVar -> ThetaPhase];
FisherCoherentPlusTher[Extremal_,ThetaPhase_, ComplexCoherent_,HilbertDim_,FieldFrequency_,Temperature_,MixConstant_] :=
Sum[((PderivCoherentPlusTher[Extremal,ThetaPhase, Iterator, ComplexCoherent,HilbertDim,FieldFrequency,Temperature,MixConstant])^2)/Tr[CoherentPlusThermalState[HilbertDim -1, ThetaPhase, ComplexCoherent, FieldFrequency, Temperature, MixConstant].Extremal[[Iterator]]], 
{Iterator, 1, Length[Extremal]}];
PderivQubit[Extremal_,ThetaAngle_, Iterator_, EtaAngle_] := Module[{DerivativeVar}, D[Tr[QBDM[EtaAngle, DerivativeVar].Extremal[[Iterator]]], DerivativeVar] /. {DerivativeVar -> ThetaAngle}];
FisherQubit[Extremal_,ThetaAngle_, EtaAngle_] := Sum[((PderivQubit[Extremal,ThetaAngle, Iterator, EtaAngle])^2)/Tr[QBDM[EtaAngle, ThetaAngle].Extremal[[Iterator]]], {Iterator, 1, Length[Extremal]}];
NormCnsnt = 1/NIntegrate[Exp[-(Theta - \[Pi])^2/(2*(\[Pi]/4)^2)], {Theta, 0, 2*\[Pi]}];
Gaussianpdf[Theta_, Sigma_, Gamma_] := NormCnsnt*Exp[-(Theta - Gamma)^2/(2*Sigma^2)];
SobraGaussian = NIntegrate[Power[D[Gaussianpdf[Theta, \[Pi]/4, \[Pi]], Theta], 2]/Gaussianpdf[Theta,\[Pi]/4, \[Pi]], {Theta, 0, 2*\[Pi]}];
(*}}}*)

End[]
EndPackage[]
