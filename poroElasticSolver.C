/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    poroElasticSolver is a modified version of solidDisplacementFoam

Description
    Transient segregated finite-volume solver of the Biot equations for
    linear-elastic, small-strain deformation of a solid skeleton coupled with 
    pore water flow and pressure governed by Darcy's law.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D and the pore water pressure p. 
    Also generating the stress tensor field sigma.
	
Author
	Johan Roenby, DHI Water & Environment

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readMaterialProperties.H"
    #include "readPoroElasticControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration: " << runTime.value() << nl << endl;

        #include "readPoroElasticControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;
        scalar pResidual = 0;

        do
        {
			volScalarField& p = pPtr();
			fvScalarMatrix pEqn
			(
				fvm::ddt(p) == fvm::laplacian(Dp, p) - fvc::div(fvc::ddt(Dp2,D))
			);

			pResidual = pEqn.solve().initialResidual();

			fvVectorMatrix DEqn
			(
				fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
			  + divSigmaExp == fvc::grad(p)
			);

			//DEqn.setComponentReference(1, 0, vector::X, 0);
			//DEqn.setComponentReference(1, 0, vector::Z, 0);

			initialResidual = DEqn.solve().initialResidual();

			if (initialResidual < pResidual){initialResidual = pResidual;}

			volTensorField gradD = fvc::grad(D);
			sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);
			divSigmaExp = fvc::div(sigmaD - (2*mu + lambda)*gradD,"div(sigmaD)");

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        #include "calculateStress.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
