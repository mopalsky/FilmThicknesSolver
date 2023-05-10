CorrectPhi
(
    U,
    phi,
    p_rgh,
    //zamiana
    //surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
    geometricZeroField(),
    pimple
);
#include "continuityErrs.H"

