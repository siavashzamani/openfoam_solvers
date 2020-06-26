Info<< "Reading transportProperties\n" << endl;
IOdictionary transportPropertiesDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar dropT
(
    transportPropertiesDict.lookup("dropT")
);

dimensionedScalar gasT
(
    transportPropertiesDict.lookup("gasT")
);

dimensionedScalar gasP
(
    transportPropertiesDict.lookup("gasP")
);

dimensionedScalar gasNu
(
    transportPropertiesDict.lookup("gasNu")
);

dimensionedScalar gasU
(
    transportPropertiesDict.lookup("gasU")
);

dimensionedScalar vaporM
(
    transportPropertiesDict.lookup("vaporM")
);

dimensionedScalar Runiversal
(
    transportPropertiesDict.lookup("Runiversal")
);

dimensionedScalar liquidNu
(
    transportPropertiesDict.lookup("liquidNu")
);

volScalarField massFlux
(
    IOobject
    (
        "massFlux",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field pSat\n" << endl;
volScalarField pSat
(
    IOobject
    (
        "pSat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh //initial value will be overwritten in calcPSatField.H
);

#include "calcPSatField.H"

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"


Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixturePtr
(
    phaseChangeTwoPhaseMixture::New(U, phi)
);

phaseChangeTwoPhaseMixture& mixture = mixturePtr();

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();




// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"





//volVectorField GRAD_ALPHA
volScalarField GRAD_ALPHA
(
    IOobject
    (
        "GRAD_ALPHA",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField alpha1_floor
(
    IOobject
    (
        "alpha1_floor",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1
);

volScalarField alpha_Area
(
    IOobject
    (
        "alpha_Area",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1
);

volScalarField alpha_volume
(
    IOobject
    (
        "alpha_volume",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1 + alpha2
);


volScalarField alpha_Press
(
    IOobject
    (
        "alpha_Press",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1
);



//-------------------------------------------------------------------------------

volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1
);

volScalarField mu_lagrangian
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1*liquidNu.value()
);
//mu_lagrangian has been initialized manually by nu of water (1e-3)


volVectorField U_lagrangian
(
    IOobject
    (
        "U_lagrangian",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*U
);


word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicCollidingCloud kinematicCloud
(
    kinematicCloudName,
    rhoInf,
    U_lagrangian,
    mu_lagrangian,
    g
);

IOobject Hheader
(
    "H",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HPtr;

if (Hheader.typeHeaderOk<volVectorField>(true))
{
    Info<< "\nReading field H\n" << endl;

    HPtr.reset(new volVectorField (Hheader, mesh));
}

IOobject HdotGradHheader
(
    "HdotGradH",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HdotGradHPtr;

if (HdotGradHheader.typeHeaderOk<volVectorField>(true))
{
    Info<< "Reading field HdotGradH" << endl;

    HdotGradHPtr.reset(new volVectorField(HdotGradHheader, mesh));
}

#include "createNonInertialFrameFields.H"



