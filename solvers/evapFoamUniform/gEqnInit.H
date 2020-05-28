    //Centering_gravity
    const dimensionedScalar gunits("gunits", dimensionSet(0,1,-2,0,0,0,0),1.0);
    const dimensionedScalar munits("gunits", dimensionSet(0,1,0,0,0,0,0),1.0);
    const dimensionedScalar sunits("gunits", dimensionSet(0,0,-1,0,0,0,0),1.0); //this unit is added to modify the central velocity implemented for lagrangian particles
    vector DropletCenter(0.0,0.0,0.0);
    scalar dropDiam(0.0);
    scalar dropDiam_old(0.0);
    scalar time_old(0.0);
    scalar time_new(0.0);

    //Spring-Damper mechanism
    vector targetCenter(0.002,0.002,0.002);
    vector x_old(targetCenter);
    vector x_new(targetCenter);
    vector centerOfDrop(targetCenter);
    vector centerOfDropVel(0.0,0.0,0.0);
    vector v_old(0.0,0.0,0.0);
    vector v_new(0.0,0.0,0.0);
    vector a_old(0.0,0.0,0.0);
    vector a_new(0.0,0.0,0.0);
    vector VELVEL(0.0,0.0,0.0);
    
    scalar k_spring(1e-3);
    scalar c_damper(1e-4);
    scalar dropMass(1);

    //PID Implementation---necessary initializations
    vector error_old(0.0,0.0,0.0);
    vector PID(0.0,0.0,0.0);
    vector Int_pid(0.0,0.0,0.0);

    //coefficients
    scalar Kp_pid(0.1);
    scalar Kd_pid(0.1);
    scalar Ki_pid(0.01);
