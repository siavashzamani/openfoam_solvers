{

volScalarField alpha1_floor_((alpha1));
float k = 0;
forAll(alpha1_floor_,celli){
	alpha1_floor_[celli] = round(alpha1_floor_[celli]);
	if(alpha1_floor_[celli]==1){k=k+1;}
}
reduce(k, sumOp<float>());
alpha1_floor = alpha1_floor_;

volScalarField pos_X(mesh.C().component(vector::X)*alpha1_floor); //weighted coordinates-x
volScalarField pos_Y(mesh.C().component(vector::Y)*alpha1_floor); //weighted coordinates-y
volScalarField pos_Z(mesh.C().component(vector::Z)*alpha1_floor); //weighted coordinates-z

vector centerOfDrop(gSum(pos_X)/k,gSum(pos_Y)/k,gSum(pos_Z)/k);
vector targetCenter(0.002,0.002,0.002);
scalar dt_(runTime.deltaT().value());
scalar cell_volume(mesh.V()[1]);
scalar distance_threshold(Foam::pow(cell_volume,1.0/3.0)*0.01);
scalar distanceFromCenter(Foam::sqrt(sqr(centerOfDrop[0]-targetCenter[0])+sqr(centerOfDrop[1]-targetCenter[1])+sqr(centerOfDrop[2]-targetCenter[2])));

DropletCenter = centerOfDrop;

if(distanceFromCenter>distance_threshold){
        g = gunits*2.0*(targetCenter-centerOfDrop)/(10.0*sqr(dt_));
        
	//Info << "g value is" << g << endl;
        //cin >> k;
}


        //g = gunits*9.81*(targetCenter-centerOfDrop)/distanceFromCenter;



//-----Lagrangian Velocity :: Virtual Velocity---------------

/*
volScalarField alpha1_abs((alpha1/(alpha1+0.01)));
forAll(alpha1_abs, celli)
{
    alpha1_abs[celli] = round(alpha1_abs[celli]);
}

volScalarField alpha2_abs(1.0 - alpha1_abs);

const volVectorField& CellCenter = mesh.C();

//volVectorField VirtualVelocity(U_units*(centerOfDrop - CellCenter)*alpha2_abs);

U_lagrangian = alpha1_abs*U; //+ VirtualVelocity;
*/



}
