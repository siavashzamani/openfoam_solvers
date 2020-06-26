{


//calculating Diffusion Coefficient: Temperature [K], Pressure [pascal]
double Diff_H2O_Air;
double ref_Patm;
ref_Patm = gasP.value()*10.0/1.01325e6;
if(gasT.value() < 450.0){
	Diff_H2O_Air = 0.0001*(1.87e-6/ref_Patm)*(Foam::pow(gasT.value(),2.072));
} else {
	Diff_H2O_Air = 0.0001*(27.5e-6/ref_Patm)*(Foam::pow(gasT.value(),1.632));
}


//Calculating Schmidt Number: Sc=gasNu/Diff_H2O_Air
double Sc_ = gasNu.value()/Diff_H2O_Air;
//Calculating Reynolds Number: Re = gasU*dropDiam/gasNu
double Re_ = gasU.value()*dropDiam/gasNu.value();
//calculating Sherwood Number
double Sh_ = 2.0 + 0.6*Foam::pow(Re_,0.5)*Foam::pow(Sc_,1.0/3.0);
//--------------------------A NOTE ON UNITS----------------------

//vaporM: vapor Molecular Weight [kg/mol]: 18*0.001 [kg/mol] water vapor ---> multiply by 0.001
//Diff_H2O: binary coefficient of diffusion [m^2/s]: 0.0001*0.3988 [m^2/s] water vapor ---> multiply by 0.0001
//Runiversal: universal gas constant [kg.m^2/(K.mol.s^2)]

//----------------------------------------------------------------
//calculating Coefficient of Mass Transfer: hm_
double hm_ = Diff_H2O_Air*Sh_/dropDiam;

//calculating Vapor Density: rho_vapor
double rho_vapor = vaporM.value()*pSat_H2O/(Runiversal.value()*dropT.value());


//calculating Mass Flux (kg.m^3/(m^2.s.kg))
massFlux = 10.*hm_*rho_vapor/rho1.value();



//to have no evaporation before the set time
if(runTime.value() < 5.0){
	massFlux = massFlux*0.0;
}





//double C_0 = 8.0*Diff_H2O_Air*vaporM.value()*pSat_H2O/(rho1.value()*Runiversal.value()*dropT.value());
//int dummy;
//std::cout.precision(18);
//std::cout << "C_0: " << std::fixed << C_0 << endl;
//cin >> dummy;


}
