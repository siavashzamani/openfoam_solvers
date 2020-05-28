{


double Coeff[8] = {-7.419242, 2.97221e-1,-1.155286e-1, 8.68563e-3, 1.094098e-3,-4.39993e-3, 2.520658e-3, -5.218684e-4};
double ref_T = 338.15;
double ref_P = 0.01;
double CP_T = 647.286;
double CP_P = 22.089e6;

//temporary variables
double sum1 = 0;
double val = ref_P*(dropT.value() - ref_T);
for(int i=0; i<8; i++){
	sum1 = sum1 + Coeff[i]*(pow(val,i));
}

val = sum1*((CP_T/dropT.value()) - 1.0);

//calculate the saturation pressure in pascals
pSat_H2O = CP_P*Foam::exp(val);




}