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

Application
    interPhaseChangeFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids with phase-change
    (e.g. cavitation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "MCCA_functions.H"
#include "basicKinematicCollidingCloud.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

		    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"


    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    //bool writeFile{true};//Added by Sia for writing custom data
    #include <iostream>
    #include <fstream>
    std::ofstream dataOutput("dataOutput.txt");
    //dataOutput.open("dataOutput.txt");
    dataOutput << "Time \t mass \t totalVolume \t DropletCenter(x) \t DropletCenter(y) \t DropletCenter(z)" << std::endl;
    dataOutput.flush();
    float t_star = 0.0;
    //float Delta_T_write;
    //int ttt;

    #include "gEqnInit.H"

    //calculate the vapor density
    scalar pSat_H2O(0.0);
    #include "pSaturation.H"

    //calculate pressures
    double P_dyn = 0.5*rho2.value()*Foam::pow(gasU.value(),2.0); //dynamic pressure
    double P_stag = 101325; //upstream pressure in pascals    



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {



		#include "gEqn.H"
		#include "gh.H"
		//g = gunits*vector(0,-1,0);
        //#include "gh.H"


		//Info << "g is:" << g << endl;
        //std::cin >> ttt;
        #include "readDyMControls.H"


        // Store divU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        volScalarField divU("divU0", fvc::div(fvc::absolute(phi, U)));

        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"
        #include "AreaEstimate.H"
        #include "massFlux.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //kinematicCloud.storeGlobalPositions();



        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {


            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;
                    Info << "inner loop->g is:" << g << endl;


                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            divU = fvc::div(fvc::absolute(phi, U));

            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            );

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            //#include "TEqn.H"
            //#include "calcPSatField.H"
            #include "writeOutput.H"


            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        //updating mu for lagrangian calculations
        #include "mu_lagrangian.H"

        kinematicCloud.evolve();


        runTime.write();
	//area estimation module

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }


    Info<< "End\n" << endl;
    //dataOutput.close();


    return 0;
}


// ************************************************************************* //
