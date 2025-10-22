/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "functionObjectTemplate.H"
#define namespaceFoam  // Suppress <using namespace Foam;>
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relVelocityFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    relVelocityFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = aef22c59d3d1023dd90795d029b22ed3b6111070
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void relVelocity_aef22c59d3d1023dd90795d029b22ed3b6111070(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::fvMesh&
Foam::relVelocityFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
relVelocityFunctionObject::
relVelocityFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
relVelocityFunctionObject::
~relVelocityFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool
Foam::
relVelocityFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        printMessage("read relVelocity");
    }

//{{{ begin code
    #line 44 "/home/cfd86/openfoam/tutorials/incompressible/pimpleFoam/RAS/rotatingFanInRoom/system/controlDict/functions/relVelocity"
const dictionary& coeffs = dict.optionalSubDict("coeffs");
        const dictionary& context = this->codeContext();

        origin = coeffs.get<vector>("origin");

        omega =
        (
            // speed
            (
                coeffs.found("rpm")
              ? degToRad(coeffs.get<scalar>("rpm") / 60.0)
              : coeffs.get<scalar>("omega")
            )
            // axis
          * normalised(coeffs.getOrDefault<vector>("axis", vector(0,0,1)))
        );

        if (!coeffs.readIfPresent("zones", zoneNames))
        {
            if (coeffs.found("cellZone"))
            {
                zoneNames.resize(1);
                coeffs.readEntry("cellZone", zoneNames[0]);
            }
        }

        if (context.getOrDefault<bool>("verbose", false))
        {
            Log<< "Relative velocity at origin " << origin << "\n";
        }
//}}} end code

    return true;
}


bool
Foam::
relVelocityFunctionObject::execute()
{
    if (false)
    {
        printMessage("execute relVelocity");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool
Foam::
relVelocityFunctionObject::write()
{
    if (false)
    {
        printMessage("write relVelocity");
    }

//{{{ begin code
    #line 78 "/home/cfd86/openfoam/tutorials/incompressible/pimpleFoam/RAS/rotatingFanInRoom/system/controlDict/functions/relVelocity"
const dictionary& context = this->codeContext();

        if (context.getOrDefault<bool>("verbose", false))
        {
            Log<< "Calculate relative velocity\n";
        }

        const auto& cc = mesh().C();
        const auto& U = mesh().lookupObject<volVectorField>("U");

        auto trelVel = volVectorField::New
        (
            "relVelocity",
            mesh(),
            dimensionedVector(dimVelocity, Zero),
            fvPatchVectorField::zeroGradientType()
        );
        auto& relVel = trelVel.ref();
        auto& relVelField = relVel.primitiveFieldRef();

        if (zoneNames.empty())
        {
            for (label celli = 0; celli < mesh().nCells(); ++celli)
            {
                relVelField[celli] = U[celli] - (omega ^ (cc[celli] - origin));
            }
        }
        else
        {
            for (const label celli : mesh().cellZones().selection(zoneNames))
            {
                relVelField[celli] = U[celli] - (omega ^ (cc[celli] - origin));
            }
        }

        relVel.correctBoundaryConditions();
        relVel.write();
//}}} end code

    return true;
}


bool
Foam::
relVelocityFunctionObject::end()
{
    if (false)
    {
        printMessage("end relVelocity");
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// ************************************************************************* //

