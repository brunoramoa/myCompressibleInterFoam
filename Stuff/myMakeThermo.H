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

InClass
    Foam::myFluidThermo

Description
    Macros for creating basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef myMakeThermo_H
#define myMakeThermo_H

#include "myFluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)          \
                                                                               \
    typedef                                                                    \
        Transport                                                              \
        <                                                                      \
            species::thermo                                                    \
            <                                                                  \
                Thermo                                                         \
                <                                                              \
                    EqnOfState                                                 \
                    <                                                          \
                        Specie                                                 \
                    >                                                          \
                >,                                                             \
                Type                                                           \
            >                                                                  \
        >                                                                      \
        Transport##Type##Thermo##EqnOfState##Specie


#define defineThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)       \
                                                                               \
    typedef                                                                    \
        CThermo                                                                \
        <                                                                      \
            BaseThermo,                                                        \
            Mixture<ThermoPhys>                                                \
        >                                                                      \
        CThermo##Mixture##ThermoPhys;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        CThermo##Mixture##ThermoPhys,                                          \
        (#CThermo"<" + Mixture<ThermoPhys>::typeName() + ">").c_str(),         \
        0                                                                      \
    )


#define addThermoPhysicsThermo(BaseThermo,CThermoMixtureThermoPhys)            \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseThermo,                                                            \
        CThermoMixtureThermoPhys,                                              \
        fvMesh                                                                 \
    );                                                                         \


#define myMakeThermoPhysicsThermo(BaseThermo,CThermo,Mixture,ThermoPhys)         \
                                                                               \
    defineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
                                                                               \
    addThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


#define myMakeThermoPhysicsThermos(BaseThermo,CThermo,Mixture,ThermoPhys)        \
                                                                               \
    defineThermoPhysicsThermo(BaseThermo, CThermo, Mixture, ThermoPhys);       \
                                                                               \
    addThermoPhysicsThermo(myBasicThermo, CThermo##Mixture##ThermoPhys);         \
    addThermoPhysicsThermo(myFluidThermo, CThermo##Mixture##ThermoPhys);         \
    addThermoPhysicsThermo(BaseThermo, CThermo##Mixture##ThermoPhys)


#define myMakeThermo(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    myMakeThermoPhysicsThermo                                                    \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )


#define myMakeThermos(BaseThermo,CThermo,Mixture,Transport,Type,Thermo,EqnOfState,Specie) \
                                                                               \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie);             \
                                                                               \
    myMakeThermoPhysicsThermos                                                   \
    (                                                                          \
        BaseThermo,                                                            \
        CThermo,                                                               \
        Mixture,                                                               \
        Transport##Type##Thermo##EqnOfState##Specie                            \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
