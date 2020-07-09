{

        //Writing data into output file
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

        //Pair<tmp<volScalarField>> mdotAlpha_ = mixture.mDotAlphal();
        //const volScalarField& mdotAlpha_ref = mdotAlpha_[1]();
        float t = runTime.value();
        const scalarField& Volumes = mesh.V();

        
        if(fabs(t-t_star) > 0.0001 || t_star ==0.0)
	{
        	dataOutput<<t<< "\t" << gSum(alpha1*Volumes)*rho1.value() << "\t"<< gSum(alpha1*Volumes) << "\t" << DropletCenter[0] << "\t" << DropletCenter[1] <<"\t" << DropletCenter[2] << std::endl;
        	dataOutput.flush();
        	
                //Info << "inside the loop:" << gSum(filter_alpha) << endl;
                //cin >> ttt;
                t_star = t;
        }
        //int tt;
        //cin >> tt;
	
	
}
