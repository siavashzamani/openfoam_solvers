/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Evaporation.H"
#include "addToRunTimeSelectionTable.H"
#include <iostream>





// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Evaporation, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Evaporation, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Evaporation::Evaporation
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi)


{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::Evaporation::mDotAlphal() const
{
    const volScalarField& alpha_Area = alpha1().db().lookupObject<volScalarField>("alpha_Area");
    const volScalarField& alpha_volume = alpha1().db().lookupObject<volScalarField>("alpha_volume");
    const volScalarField& massFlux = alpha1().db().lookupObject<volScalarField>("massFlux"); //kg/(m^2*s)
    const scalar rhoGas(rho2().value()); //(kg/m^3)
    //const dimensionedScalar pdummy2 ("pdummy2", dimensionSet(1,-3,-1,0,0,0,0), -0.05);
    //const dimensionedScalar U_("U_", dimensionSet(0,0,-1,0,0,0,0), 1.5e-5);
    //const dimensionedScalar rho_gas("rho_gas", dimensionSet(1,-3,0,0,0,0,0), 1.0);
    const volScalarField alphaSupreme(alpha1()+alpha2());
    const volScalarField alpha1_dividable(alpha1() + SMALL);
    const dimensionedScalar pdummy ("pdummy", dimensionSet(1,-3,-1,0,0,0,0), 0.0);
    const dimensionedScalar dimdummy("dimdummy", dimensionSet(1,-3,-1,0,0,0,0), 1.0);
    
    
    return Pair<tmp<volScalarField>>
    (
        alphaSupreme*pdummy,
        //alphaSupreme*pdummy2
        -1.0*dimdummy*rhoGas*massFlux*alpha_Area/((alpha_volume+SMALL)*(alpha1_dividable))
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixtures::Evaporation::mDotP() const
{
    const volScalarField alphaSupreme(alpha1()+alpha2());
    const dimensionedScalar pdummy1 ("pdummy1", dimensionSet(0,-2,1,0,0,0,0), 0.0); 


    return Pair<tmp<volScalarField>>
    (
        alphaSupreme*pdummy1,
        alphaSupreme*pdummy1
    );
}


void Foam::phaseChangeTwoPhaseMixtures::Evaporation::correct()
{
    phaseChangeTwoPhaseMixture::correct();
}



bool Foam::phaseChangeTwoPhaseMixtures::Evaporation::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        //phaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}




// ************************************************************************* //
