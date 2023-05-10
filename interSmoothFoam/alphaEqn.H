{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    tmp<fv::ddtScheme<scalar>> ddtAlpha(
        fv::ddtScheme<scalar>::New(
            mesh,
            mesh.ddtScheme("ddt(alpha)")));

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    if (
        isType<fv::EulerDdtScheme<scalar>>(ddtAlpha()) || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha()))
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha()))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
                << "Sub-cycling is not supported "
                   "with the CrankNicolson ddt scheme"
                << exit(FatalError);
        }

        ocCoeff =
            refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha())
                .ocCoeff();
    }
    else
    {
        FatalErrorInFunction
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    scalar cnCoeff = 1.0 / (1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha() * mag(phi / mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha() * icAlpha) * fvc::interpolate(mag(U));
    }

    surfaceScalarField::Boundary &phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField &phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff * phi + (1.0 - cnCoeff) * phi.oldTime();
    }


    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
       // - fvm::Sp(fvc::ddt(dimensionedScalar("1", dimless, 1), mesh)
       //           + fvc::div(phiCN), alpha1)
         ==
            Su + fvm::Sp(Sp + divU, alpha1)
        );


        Info << "Phase-1 volume fraction = "
             << alpha1.weightedAverage(mesh.Vsc()).value()
             << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
             << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
             << endl;

        tmp<surfaceScalarField>  talphaPhiUD(alpha1Eqn.flux());
        alphaPhi = talphaPhiUD();

        // Nie wiem czy to jest dobrze
      if (alphaApplyPrevCorr &&  talphaPhiCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                alpha1,
                alphaPhi,
                 talphaPhiCorr0.ref(),
                oneField(),
                zeroField()
            );

            alphaPhi +=  talphaPhiCorr0();
        }

        // Cache the upwind-flux
         talphaPhiCorr0 = talphaPhiUD;

        alpha2 = 1.0 - alpha1;
        interface.correct();
    }

    volScalarField alpha10("alpha10", alpha1);

    for (int aCorr = 0; aCorr < nAlphaCorr; aCorr++)
    {

        #include "alphaSuSp.H"

        surfaceScalarField phir(phic * interface.nHatf());

        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                cnCoeff*alpha1 + (1.0 - cnCoeff)*alpha1.oldTime(),
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(alphaPhiUn - alphaPhi);
            // volScalarField alpha10("alpha10", alpha1);
            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            // MULES::correct(alpha1, alphaPhiUn, talphaPhiCorr.ref(), 1, 0);
            MULES::correct(
                geometricOneField(),
                alpha1,
                alphaPhiUn,
                talphaPhiCorr.ref(),
                oneField(),
                zeroField());

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5 * alpha1 + 0.5 * alpha10;
                alphaPhi += 0.5 * talphaPhiCorr();
            }
        }
        else
        {
            alphaPhi = alphaPhiUn;
            MULES::explicitSolve(
                geometricOneField(),
                alpha1,
                phiCN,
                alphaPhi,
                Sp,
                (Su + divU*min(alpha1(), scalar(1)))(),
                oneField(),
                zeroField()
                );
                
        }

        alpha2 = 1.0 - alpha1;

        // nie wiem czy to dobrze jest
        // brane jest z interFoam
        interface.correct();
    }


    if (alphaApplyPrevCorr && MULESCorr)
    {
         talphaPhiCorr0 = alphaPhi -  talphaPhiCorr0;
         talphaPhiCorr0.ref().rename("alphaPhi1Corr0");
    }
    else
    {
         talphaPhiCorr0.clear();
    }

    #include "rhofs.H"

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
     || word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::localEulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi*(rho1f - rho2f) + phiCN*rho2f;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi =
                (alphaPhi - (1.0 - cnCoeff)*alphaPhi.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi*(rho1f - rho2f) + phi*rho2f;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}