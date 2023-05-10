    MRF.correctBoundaryVelocity(U);




    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );
  
    /*
       fvVectorMatrix UEqn
    (
	    rho*
        (
            fvm::ddt(U) 
		  + fvm::div(phi, U)
          + MRF.DDt(U)
		)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
    );
  */


    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    //zmieniono
                    //mixture.surfaceTensionForce()
                    interface.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
