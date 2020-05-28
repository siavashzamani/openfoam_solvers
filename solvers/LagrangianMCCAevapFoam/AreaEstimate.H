{


int dummyPut;

volScalarField filter_alpha((alpha1*alpha2/(alpha1*alpha2+0.01)));
forAll(filter_alpha, celli)
{
    filter_alpha[celli] = round(filter_alpha[celli]);
}
   
GRAD_ALPHA = filter_alpha;

//GRAD_ALPHA = filter_alpha*fvc::reconstruct(fvc::snGrad(alpha1)*mesh.magSf());

//creating the normal vector to alpha interface.
volVectorField alpha_normal(fvc::reconstruct(fvc::snGrad(alpha1)));

forAll(alpha_normal, cellID)
{
	if(mag(alpha_normal[cellID]) !=0.0){
		alpha_normal[cellID] = alpha_normal[cellID]/mag(alpha_normal[cellID]);
	}

} 


const dimensionedScalar dum("dum", dimensionSet(0,4,-1,0,0,0,0),1.0);
GRAD_ALPHA1 = alpha_normal*dum*GRAD_ALPHA;


const scalarField& Volumes = mesh.V(); //total volume of drop
const scalar L_ = Foam::pow((Volumes[1]),(1.0/3.0)); //length of a single cell - assuming uniform mesh


scalar surfaceArea_total;
if(mesh.nGeometricD()==2){
//2D Averaged Surface Area Method
surfaceArea_total = 2.0*Foam::pow(constant::mathematical::pi*L_*gSum(alpha1*Volumes),0.5);
}else{
//3D Averaged Surface Area Method
surfaceArea_total = 4.0*constant::mathematical::pi*Foam::pow(3.0*gSum(alpha1*Volumes/(4.0*constant::mathematical::pi)),2.0/3.0);
}

//reseting alpha_volume value to 1
alpha_volume = alpha1+alpha2;


volScalarField Area(alpha1*0.0);
volScalarField volume_holder(alpha1+alpha2);
forAll(Area,i)
{
	if(GRAD_ALPHA[i] != 0)
	{
		double alpha_ = alpha1[i];
		//if(alpha_ > 0.5){alpha_ = 1.0 - alpha_;}
		const scalar x_ = alpha_normal[i].component(0);
		const scalar y_ = alpha_normal[i].component(1);
		const scalar z_ = alpha_normal[i].component(2);
		const scalar V_cell = mesh.V()[i];
		//const scalar L_ = Foam::pow((V_cell),(1.0/3.0));
		//const scalar theta_max = Foam::atan(L_/(2.0*alpha_*V_cell));
        

                    

		//double theta = constant::mathematical::pi/2.0 - Foam::acos(fabs(x_));
		//Weight = f(n,sqrt(alpha/0.5))
		//const scalar weight_ = 1.0 - (1.0 - Foam::sqrt(alpha_/0.5))*Foam::sin(2.0*theta);

		//Weight = f(sqrt(alpha/0.5))
		//const scalar weight_ = Foam::pow(alpha_/0.5,0.5);
		
		//Weight = f(alpha_/0.5)
		//const scalar weight_ = alpha_/0.5;


                    double theta = Foam::acos(fabs(x_/Foam::sqrt(x_*x_ + y_*y_)));
                    double phi = Foam::acos(fabs(z_));
                    double AA;
                    MCCA_Iterative_Area(alpha_, theta, phi, AA);
                    
		Area[i] = AA; //L_*L_*AA; // for exact measurement L*L must be multiplied by the calculated area for unit cube
		//Area[i] = weight_;
                    volume_holder[i] = alpha1[i]*V_cell;
                    




		
		/*
		//2D Version
		//double theta = constant::mathematical::pi/2.0 - Foam::atan(fabs(y_/x_));
		//3D Version
		double theta = constant::mathematical::pi/2.0 - Foam::asin(fabs(y_));
		const scalar Ax_ = L_*L_;
		if(theta > constant::mathematical::pi/4.0){theta = constant::mathematical::pi/2.0 - theta;}
		const scalar Area_template = Ax_*Foam::sqrt(1.0 + 1.0/Foam::pow(Foam::atan(constant::mathematical::pi/2.0 - theta),2.0));
		const scalar w_ = 1.0 - (1.0 - Foam::sqrt(alpha_/0.5)*Foam::sin(2.0*theta));
		Area[i] = Area_template*w_;
		*/



	}
}

const scalar weight_total = gSum(Area);






if(weight_total != 0.0)
{
	forAll(Area,i)
	{
		Area[i] = Area[i]*surfaceArea_total/weight_total;


	}
}


alpha_volume = volume_holder;
alpha_Area = Area;


forAll(alpha_volume,i)
{
	if(alpha_volume[i] == 0.0)
	{
		Info << "zero is encountered" << endl;
		cin >> dummyPut;
	}
}


}
