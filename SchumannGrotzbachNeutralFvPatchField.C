/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "SchumannGrotzbachNeutralFvPatchField.H"
//#include "singlePhaseTransportModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SchumannGrotzbachNeutralFvPatchField::
SchumannGrotzbachNeutralFvPatchField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(p, iF),
    kappa_(0.40),
    z0_(p.size(), 0.01),
    betaM_(15.0),
    gammaM_(4.7),
    averageType_("local")
{}


SchumannGrotzbachNeutralFvPatchField::
SchumannGrotzbachNeutralFvPatchField
(
    const SchumannGrotzbachNeutralFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    z0_(ptf.z0_, mapper),
    betaM_(ptf.betaM_),
    gammaM_(ptf.gammaM_),
    averageType_(ptf.averageType_)
{}


SchumannGrotzbachNeutralFvPatchField::
SchumannGrotzbachNeutralFvPatchField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(p, iF, dict),
    kappa_(readScalar(dict.lookup("kappa"))),
    z0_("z0", dict, p.size()),
    betaM_(readScalar(dict.lookup("betaM"))),
    gammaM_(readScalar(dict.lookup("gammaM"))),
    averageType_(dict.lookupOrDefault<word>("averageType","local"))
{}


SchumannGrotzbachNeutralFvPatchField::
SchumannGrotzbachNeutralFvPatchField
(
    const SchumannGrotzbachNeutralFvPatchField& wfpf
)
:
    fixedValueFvPatchSymmTensorField(wfpf),
    kappa_(wfpf.kappa_),
    z0_(wfpf.z0_),
    betaM_(wfpf.betaM_),
    gammaM_(wfpf.gammaM_),
    averageType_(wfpf.averageType_)
{}


SchumannGrotzbachNeutralFvPatchField::
SchumannGrotzbachNeutralFvPatchField
(
    const SchumannGrotzbachNeutralFvPatchField& wfpf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(wfpf, iF),
    kappa_(wfpf.kappa_),
    z0_(wfpf.z0_),
    betaM_(wfpf.betaM_),
    gammaM_(wfpf.gammaM_),
    averageType_(wfpf.averageType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SchumannGrotzbachNeutralFvPatchField::evaluate
(
    const Pstream::commsTypes
)
{
    // ---Get preliminary information
    //    Get face normal vectors
    const vectorField normal = patch().nf();

    //    Get face areas (individual and global sum) (note that "gSum" is used as opposed
    //    to "sum" because "gSum" is parallel-aware --- it gathers sums from each processor
    //    to which this patch belongs and makes the global sum)
    const scalarField area = patch().magSf();
    const scalar areaTotal = gSum(area);

    //    Get perpendicular distance from cell center to boundary
    const scalarField z1 = 1.0/patch().deltaCoeffs();

    vectorField loc = patch().Cf();


    // ---Get the velocity parallel to the boundary using,
    //
    //     U_|| = U_1 - ((U_1 dot n_f) * n_f), where U_||,
    //
    //    where U_|| is the parallel velocity vector, U_1 is the cell center velocity, and
    //    n_f is the surface face normal unit vector.
    //    Get the velocity in the cells adjacent to the boundary
    const fvPatchVectorField& UPatch = patch().lookupPatchField<volVectorField, vector>("U");
    vectorField UParallel = UPatch.patchInternalField();
    UParallel = UParallel - ((UParallel & normal) * normal);
    vector UParallelMean = gSum(UParallel * area) / areaTotal;
    scalar UParallelMeanMag = mag(UParallelMean);

    // ---Get the velocity parallel to the boundary in terrain-local coordinates
    vectorField UParallelP = UParallel;
    forAll(UParallelP, facei)
    {
        // Transform from the Cartesian (x,y,z) system into the local (x',y',z')
        // system

        // z' is equal to the surface normal pointing inward (negative because
        // OpenFOAM normal is outward)
        vector zP;
        zP = -normal[facei];

        // x' is pointed in the direction of the parallel resolved velocity at this cell
        vector xP;
        xP = UParallel[facei];

        // y' is orthogonal to x' and z', so it can be found with cross product
        vector yP;
        yP = zP ^ xP;

        // Transform the velocity from Cartesian to local
        UParallelP[facei] = transformVectorCartToLocal(UParallel[facei], xP, yP, zP);

      //Info << "facei = " << facei << tab
      //     << "UParallel  = " << UParallel[facei] << tab
      //     << "UParallelP = " << UParallelP << tab;

    }

    //    Get magnitudes and means of the terrain-local velocity
    scalarField UParallelPMag = mag(UParallelP);
    //scalar UParallelPMagMean = gSum(UParallelPMag * area) / areaTotal;
    vector UParallelPMean = gSum(UParallelP * area) / areaTotal;
    scalar UParallelPMeanMag = mag(UParallelPMean);


    Info << "UParallelMeanMag = " << UParallelMeanMag << tab
         << "UParallelPMeanMag = " << UParallelPMeanMag << endl;


    //    Define friction velocity
    scalar uStarMean(0.0);


    // ---Specify surface shear stresses (Schumann formulation)
    symmTensorField& Rw = *this;
    scalarField RwMag(patch().size(),0.0);

    // Form the wall shear stress tensor in terrain-local coordinates
    symmTensor RwP(symmTensor::zero);
    vector xP(vector::zero);
    vector yP(vector::zero);
    vector zP(vector::zero);


    // Info << "planarAverage" << endl;
    if(averageType_ == "planarAverage")
    {
        scalar z1Mean = gSum(z1 * area)/areaTotal;

        //    Get the average surface roughness height
        scalar z0Mean = gSum(z0_ * area)/areaTotal;

        uStarMean = (kappa_ * UParallelMeanMag) / Foam::log(z1Mean / z0Mean);

        forAll(Rw, facei)
        {
            Rw[facei].xx() = 0.0;
            Rw[facei].xy() = 0.0;
            Rw[facei].xz() = -sqr(uStarMean) * (UParallel[facei].x() / max(UParallelMeanMag, 1.0E-5));
            Rw[facei].yy() = 0.0;
            Rw[facei].yz() = -sqr(uStarMean) * (UParallel[facei].y() / max(UParallelMeanMag, 1.0E-5));
            Rw[facei].zz() = 0.0;

            // Get the magnitude of the surface stress vector
            RwMag[facei] = Foam::sqrt(Foam::sqr(Rw[facei].xz()) + Foam::sqr(Rw[facei].yz()));
        }
    }

    //          Info << "local" << endl;
    else if (averageType_ == "local")
    {
        scalarField uStar(patch().size(),0.0);

        forAll(uStar, facei)
        {
            uStar[facei] = (kappa_ * UParallelPMag[facei]) / Foam::log(z1[facei] / z0_[facei]);

            RwP.xx() = 0.0;
            RwP.xy() = 0.0;
            RwP.xz() = -sqr(uStar[facei]) * (UParallelP[facei].x() / max(UParallelPMag[facei], 1.0E-5));
            RwP.yy() = 0.0;
            RwP.yz() = -sqr(uStar[facei]) * (UParallelP[facei].y() / max(UParallelPMag[facei], 1.0E-5));
            RwP.zz() = 0.0;


          // Transform from the terrain-local into the Cartesian system

          // z' is equal to the surface normal pointing inward (negative because
          // OpenFOAM normal is outward)
          zP = -normal[facei];

          // x' is pointed in the direction of the parallel resolved velocity at this cell
          xP = UParallel[facei];

          // y' is orthogonal to x' and z', so it can be found with cross product
          yP = zP ^ xP;

          // Perform the transformation
          Rw[facei] = transformSymmTensorLocalToCart(RwP, xP, yP, zP);

          // Get the magnitude of the surface stress vector
          RwMag[facei] = Foam::sqrt(Foam::sqr(RwP.xz()) + Foam::sqr(RwP.yz()));
        }

        //    Get average friction velocity
        uStarMean = gSum(uStar * area) / areaTotal;
    }


    symmTensor RwMean = gSum(Rw * area) / areaTotal;
    scalar RwMeanMag = Foam::sqrt(Foam::sqr(RwMean.xz()) + Foam::sqr(RwMean.yz()));
    scalar RwMagMean = gSum(RwMag * area) / areaTotal;

    Info << "RwMagMean = " << RwMagMean << tab
         << "RwMeanMag = " << RwMeanMag << tab
         << "sqrt(RwMagMean) = " << Foam::sqrt(RwMagMean) << tab
         << "sqrt(RwMeanMag) = " << Foam::sqrt(RwMeanMag) << tab
         << "uStarMean = " << uStarMean << endl;
}

vector SchumannGrotzbachNeutralFvPatchField::transformVectorCartToLocal
(
    vector v,
    vector xP,
    vector yP,
    vector zP
)
{
    // Transform from the Cartesian (x,y,z) system into the local (x',y',z')
    // system
    //
    //    x' is aligned with the flow
    //    y' is the cross product of z' and x'
    //    z' is in the boundary face normal direction
    //
    // These vectors are unit vectors.  The vectors make up the rows of
    // the rotation matrix T', which rotates from (x,y,z) to (x',y',z').

    // z' is equal to the surface normal pointing inward (negative because
    // OpenFOAM normal is outward)
    scalar zPMag;
    zPMag = mag(zP);
    if (zPMag != 0.0)
    {
       zP = zP/zPMag;
    }

    // x' is pointed in the direction of the parallel resolved velocity at this cell
    scalar xPMag;
    xPMag = mag(xP);
    if (xPMag != 0.0)
    {
       xP = xP/xPMag;
    }

    // y' is orthogonal to x' and z', so it can be found with cross product
    scalar yPMag;
    yPMag = mag(yP);
    if (yPMag != 0.0)
    {
       yP = yP/yPMag;
    }

    // Create T'
    tensor TP;
    TP.xx() = xP.x();
    TP.xy() = xP.y();
    TP.xz() = xP.z();
    TP.yx() = yP.x();
    TP.yy() = yP.y();
    TP.yz() = yP.z();
    TP.zx() = zP.x();
    TP.zy() = zP.y();
    TP.zz() = zP.z();

    // Transform the vector from Cartesian to local
    vector vP = TP & v;

    //Info << "xP = " << xP << tab
    //     << "yP = " << yP << tab
    //     << "zP = " << zP << tab
    //     << "v  = " << v  << tab
    //     << "vP = " << vP << endl;

    return vP;
}

symmTensor SchumannGrotzbachNeutralFvPatchField::transformSymmTensorLocalToCart
(
    symmTensor SP,
    vector xP,
    vector yP,
    vector zP
)
{
    // Transform from the local (x',y',z') system into the Cartesian (x,y,z)
    // system
    //
    //    x' is aligned with the flow
    //    y' is the cross product of z' and x'
    //    z' is in the boundary face normal direction
    //
    // These vectors are unit vectors.  The vectors make up the rows of
    // the rotation matrix, T', which rotates from (x,y,z) to (x',y',z').
    // T'^-1 = T transforms from (x',y',z') to (x,y,z).  Note that T' is
    // such that the (T')^-1 = transpose(T') because it is made up of
    // orthogonal basis vectors.  A tensor can be transformed from (x',y',z')
    // to (x,y,z) using S_ij = T * S_ij' * T'.

    // z' is equal to the surface normal pointing inward (negative because
    // OpenFOAM normal is outward)
    scalar zPMag;
    zPMag = mag(zP);
    if (zPMag != 0.0)
    {
       zP = zP/zPMag;
    }

    // x' is pointed in the direction of the parallel resolved velocity at this cell
    scalar xPMag;
    xPMag = mag(xP);
    if (xPMag != 0.0)
    {
       xP = xP/xPMag;
    }

    // y' is orthogonal to x' and z', so it can be found with cross product
    scalar yPMag;
    yPMag = mag(yP);
    if (yPMag != 0.0)
    {
       yP = yP/yPMag;
    }

    // Create T'
    tensor TP;
    TP.xx() = xP.x();
    TP.xy() = xP.y();
    TP.xz() = xP.z();
    TP.yx() = yP.x();
    TP.yy() = yP.y();
    TP.yz() = yP.z();
    TP.zx() = zP.x();
    TP.zy() = zP.y();
    TP.zz() = zP.z();

    // Create T by transposing T' to get T'^-1
    tensor T;
    T = TP.T();

    // Transform the symmetric tensor by doing S_ij = T * S_ij' * T'
    tensor Ss;
    Ss = T & SP & TP;

    // OpenFOAM will not allow a tensor & tensor operation resulting in a symmTensor,
    // so the work around was to make Ss a tensor, even though it ends up being a
    // symmetric one, and then copy the relevant components into symmTensor S
    symmTensor S(symmTensor::zero);
    S.xx() = Ss.xx();
    S.xy() = Ss.xy();
    S.xz() = Ss.xz();
    S.yy() = Ss.yy();
    S.yz() = Ss.yz();
    S.zz() = Ss.zz();

    //Info << "SP = " << SP << endl
    //     << "S  = " << S  << endl;

    return S;
}



void SchumannGrotzbachNeutralFvPatchField::write(Ostream& os) const
{
    fvPatchField<symmTensor>::write(os);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    z0_.writeEntry("z0", os);
    os.writeKeyword("betaM") << betaM_ << token::END_STATEMENT << nl;
    os.writeKeyword("gammaM") << gammaM_ << token::END_STATEMENT << nl;
    os.writeKeyword("averageType") << averageType_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchSymmTensorField,
    SchumannGrotzbachNeutralFvPatchField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
