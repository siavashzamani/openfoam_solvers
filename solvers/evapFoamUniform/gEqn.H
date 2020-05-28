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

volScalarField U_X(U.component(0)*alpha1_floor);
volScalarField U_Y(U.component(1)*alpha1_floor);
volScalarField U_Z(U.component(2)*alpha1_floor);


x_old = centerOfDrop;
v_old = centerOfDropVel;
centerOfDrop[0] = gSum(pos_X)/k;
centerOfDrop[1] = gSum(pos_Y)/k;
centerOfDrop[2] = gSum(pos_Z)/k;

centerOfDropVel[0] = gSum(U_X)/k;
centerOfDropVel[1] = gSum(U_Y)/k;
centerOfDropVel[2] = gSum(U_Z)/k;

x_new = centerOfDrop;
v_new = centerOfDropVel;


//PID implementation----------------------------------------------
scalar dt_(runTime.deltaT().value());
vector error_pid(targetCenter - x_new);
vector p_pid(error_pid);
vector d_pid((error_pid - error_old)/dt_);
vector i_pid((error_pid + error_old)*dt_/2.0);
Int_pid = Int_pid + i_pid;

PID = Kp_pid*p_pid + Kd_pid*d_pid + Ki_pid*Int_pid;


//PID implementation----------------------------------------------


scalar cell_volume(mesh.V()[1]);
scalar distance_threshold(Foam::pow(cell_volume,1.0/3.0)*0.01);
scalar distanceFromCenter(Foam::sqrt(sqr(centerOfDrop[0]-targetCenter[0])+sqr(centerOfDrop[1]-targetCenter[1])+sqr(centerOfDrop[2]-targetCenter[2])));

DropletCenter = centerOfDrop;

vector springX(targetCenter- centerOfDrop);
//vector springV((x_new-x_old)/(dt_));
vector springV(centerOfDropVel);
vector springA((v_new-v_old)/(dt_));
const scalarField& Volumes = mesh.V();
dropMass = gSum(Volumes*alpha1)*rho1.value();
VELVEL = springV;


vector gravity_((k_spring/dropMass)*springX - (c_damper/dropMass)*springV + PID/dropMass);

g = gunits*gravity_;



error_old = error_pid;

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
