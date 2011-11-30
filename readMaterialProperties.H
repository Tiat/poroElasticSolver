    Info<< "Reading material properties\n" << endl;

    IOdictionary materialProperties
    (
        IOobject
        (
            "materialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar E(materialProperties.lookup("E"));
    dimensionedScalar nu(materialProperties.lookup("nu"));
    dimensionedScalar k(materialProperties.lookup("k"));
    dimensionedScalar K(materialProperties.lookup("K"));
    dimensionedScalar Sr(materialProperties.lookup("Sr"));
    dimensionedScalar p0(materialProperties.lookup("p0"));
    dimensionedScalar gamma(materialProperties.lookup("gamma"));
    dimensionedScalar n(materialProperties.lookup("n"));

    dimensionedScalar Kprime = 1.0/(1.0/K + (1.0-Sr)/p0);

    Info<< "Calculating Lame's coefficients\n" << endl;

    dimensionedScalar mu = E/(2.0*(1.0 + nu));
    dimensionedScalar lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));

    Switch planeStress(materialProperties.lookup("planeStress"));

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }

    Info<< "mu = " << mu.value() << " Pa\n";
    Info<< "lambda = " << lambda.value() << " Pa\n";

    dimensionedScalar Dp
    (
        "Dp",
        dimensionSet(0, 2, -1 , 0, 0),
        (k/gamma*Kprime/n).value()
    );

    dimensionedScalar Dp2
    (
        "Dp2",
        dimensionSet(1, -1, -2 , 0, 0),
        (Kprime/n).value()
    );

    Info<< "Dp = " << Dp.value() << " m^2/s \n";
    Info<< "Dp2 = " << Dp2.value() << " kg/m/s^2\n";

