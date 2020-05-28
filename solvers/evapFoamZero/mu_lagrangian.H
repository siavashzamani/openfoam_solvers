{


volScalarField alpha1_abs((alpha1/(alpha1+0.01)));
forAll(alpha1_abs, celli)
{
    alpha1_abs[celli] = round(alpha1_abs[celli]);
}

volScalarField alpha2_abs(1.0 - alpha1_abs);

rhoInf = alpha1_abs*rho1 + alpha2_abs*rho2;

//volVectorField U_alpha(U*alpha1_abs);


volVectorField VirtualVelocity((munits*DropletCenter - mesh.C())*alpha2_abs);
time_old = time_new;
time_new = runTime.value();
scalar virtualVelocityMagnitude = fabs(dropDiam_old - dropDiam)/(time_new - time_old); //mag(gMax(U_alpha));

//normalize virtual velocity -- at this point, you can multiply the normalized velocity by an appropriate magnitude to match surface tension conditions.
forAll(VirtualVelocity, celli){
	VirtualVelocity[celli] = virtualVelocityMagnitude*VirtualVelocity[celli]/(mag(VirtualVelocity[celli])+VSMALL);
}

U_lagrangian = U; //alpha1_abs*U + sunits*VirtualVelocity;
//this velocity should be scaled based on the surface tension force acting on solid particles.
mu_lagrangian = liquidNu.value()*alpha1_abs*rho1 + gasNu.value()*alpha2_abs*rho2;


}
