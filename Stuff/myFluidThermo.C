/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "myFluidThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(myFluidThermo, 0);
    defineRunTimeSelectionTable(myFluidThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myFluidThermo::myFluidThermo(const fvMesh& mesh, const word& phaseName)
:
    myBasicThermo(mesh, phaseName)
{}



Foam::myFluidThermo::myFluidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    myBasicThermo(mesh, dict, phaseName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::myFluidThermo> Foam::myFluidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return myBasicThermo::New<myFluidThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myFluidThermo::~myFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::myFluidThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::myFluidThermo::nu(const label patchi) const
{
    return mu(patchi)/rho(patchi);
}


// ************************************************************************* //
