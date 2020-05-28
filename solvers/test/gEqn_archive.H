{
	 
	 /*
	 volScalarField alpha_p(alpha1-max(alpha1));
	 volScalarField alpha_p1(alpha_p/(alpha_p+1.0e-200));
         volScalarField alpha1_final(alpha1*alpha_p1);
	 volScalarField alpha1_normal(cmptDivide(alpha1_final,(alpha1_final+1.0e-20)));
	 GRAD_ALPHA = alpha1*alpha2/(alpha1*alpha2+0.01);
	 */
/*
	 volScalarField filter_alpha((alpha1*alpha2/(alpha1*alpha2+0.01)));
	 forAll(filter_alpha, celli){
	 filter_alpha[celli] = round(filter_alpha[celli]);
	 }


         if(gSum(filter_alpha)<0.0){
		volScalarField alpha1dummy(alpha1);
		volScalarField alpha2dummy(alpha2);
		forAll(alpha1dummy,celli){
			if(alpha1dummy[celli]>1.0){
				alpha1dummy[celli]=1.0;
				alpha2dummy[celli]=0.0;
			}
		}
		filter_alpha = (alpha1dummy*alpha2dummy/(alpha1dummy*alpha2dummy+0.01));
		forAll(filter_alpha, celli){
			filter_alpha[celli] = round(filter_alpha[celli]);
		}
		//Info << "inside the inner loop" << gSum(filter_alpha) << endl;
                //cin >> ttt;

			
	 }

         volScalarField LAlpha1(min(max(alpha1, scalar(0)),scalar(1)));
	 volScalarField filter_alpha(LAlpha1);
	 forAll(filter_alpha,celli){
		if(filter_alpha[celli] > 0 && filter_alpha[celli] < 1){
			filter_alpha[celli] = 1.0;
		} else {
			filter_alpha[celli] = 0.0;
		}
	 }
*/        
         //Pair<tmp<volScalarField>> vDotAlphal = mixture.vDotAlphal();
         //const volScalarField& vDotcAlphal = vDotAlphal[0]();
         //const volScalarField& vDotvAlphal = vDotAlphal[1]();
         //GRAD_ALPHA = filter_alpha*fvc::reconstruct(fvc::snGrad(alpha1)*mesh.magSf()); //gradient-> volVectorField
         //GRAD_ALPHA = filter_alpha*vDotvAlphal; //this is the working line
	 //Pair<tmp<volScalarField>> mtestAlphal = mixture.mDotAlphal();
         //const volScalarField& testm = mtestAlphal[1]();
	 //volScalarField testm1(min(max(alpha1, scalar(0)),scalar(1)));
	 
         volScalarField LAlpha1(min(max(alpha1, scalar(0)),scalar(1)));
         GRAD_ALPHA = LAlpha1;

         const scalarField& Volumes = mesh.V();
         //const scalar t = mesh.time().value();
         const scalar t = runTime.value();
         const scalar dt = runTime.deltaTValue();

         Info << "sum of all alpha1:" << gSum(alpha1*Volumes) << endl;
         Info << "time is:" << t << endl;
         Info << "time step is:" << dt << endl;
         


	 /*
	 forAll(GRAD_ALPHA, celli){
		GRAD_ALPHA[celli].component(0) = round(GRAD_ALPHA[celli].component(0));
		GRAD_ALPHA[celli].component(1) = round(GRAD_ALPHA[celli].component(1));
		GRAD_ALPHA[celli].component(2) = round(GRAD_ALPHA[celli].component(2));
	 }
	 */

}
