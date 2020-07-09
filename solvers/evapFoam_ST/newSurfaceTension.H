
//surface tension as a function of distance from droplet center with dimensionSet(1,0,-2,0,0,0,0)
dimensionedScalar div_ ("div_", dimensionSet(0,1,0,0,0,0,0), 1.0e-6);
dimensionedVector DropletCenter_(dimLength, DropletCenter);
volScalarField k_STmod(ST_modUnits*sigma_STmod/(mag(mesh.C() - DropletCenter_)+div_));



//constant surface tension for the entire interface with dimensionSet(1,-1,-2,0,0,0,0)
//volScalarField k_STmod(pow(alpha1,0.0)*ST_modUnits*sigma_STmod/(dropDiam/2.0));




