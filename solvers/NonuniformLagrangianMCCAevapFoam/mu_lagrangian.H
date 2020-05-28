{


volScalarField alpha1_abs((alpha1/(alpha1+0.01)));
forAll(alpha1_abs, celli)
{
    alpha1_abs[celli] = round(alpha1_abs[celli]);
}

volScalarField alpha2_abs(1.0 - alpha1_abs);

rhoInf = alpha1_abs*rho1 + alpha2_abs*rho2*100000;

//int dummy;
//Info << "check to see this: " << (munits*DropletCenter - mesh.C()) << endl;
//std::cin >> dummy;

volVectorField VirtualVelocity((munits*DropletCenter - mesh.C())*alpha2_abs);

//normalize virtual velocity -- at this point, you can multiply the normalized velocity by an appropriate magnitude to match surface tension conditions.
forAll(VirtualVelocity, celli){
	VirtualVelocity[celli] = (1e-4)*VirtualVelocity[celli]/(mag(VirtualVelocity[celli])+VSMALL);
}

U_lagrangian = alpha1_abs*U + sunits*VirtualVelocity;
//this velocity should be scaled based on the surface tension force acting on solid particles.




//U_lagrangian = alpha1_abs*U;
//U_lagrangian = alpha1*U;



mu_lagrangian = 1e-3*alpha1_abs*rho1 + alpha2_abs*rho2*1.48e-05;


}
