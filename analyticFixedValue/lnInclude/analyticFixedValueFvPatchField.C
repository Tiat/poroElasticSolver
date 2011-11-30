/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "analyticFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void analyticFixedValueFvPatchField<Type>::setCurField()
{
	const vectorField x(this->patch().Cf());
	scalar t = this->db().time().value();
        this->operator==(unitField_*amplitude_*cos((waveVector_ & x) - frequency_*t));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
analyticFixedValueFvPatchField<Type>::analyticFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    unitField_(p.size()),
    amplitude_(0.0),
    frequency_(0.0),
    waveVector_(vector(1,0,0)),
    curTimeIndex_(-1)
{}


template<class Type>
analyticFixedValueFvPatchField<Type>::analyticFixedValueFvPatchField
(
    const analyticFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    unitField_(ptf.unitField_, mapper),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    waveVector_(ptf.waveVector_),
    curTimeIndex_(-1)
{}


template<class Type>
analyticFixedValueFvPatchField<Type>::analyticFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    unitField_("unitField", dict, p.size()),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    waveVector_(dict.lookup("waveVector")),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(unitField_);
    }
}


template<class Type>
analyticFixedValueFvPatchField<Type>::analyticFixedValueFvPatchField
(
    const analyticFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    unitField_(ptf.unitField_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    waveVector_(ptf.waveVector_),
    curTimeIndex_(-1)
{}


template<class Type>
analyticFixedValueFvPatchField<Type>::analyticFixedValueFvPatchField
(
    const analyticFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    unitField_(ptf.unitField_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    waveVector_(ptf.waveVector_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void analyticFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    unitField_.autoMap(m);
}


template<class Type>
void analyticFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const analyticFixedValueFvPatchField<Type>& tiptf =
        refCast<const analyticFixedValueFvPatchField<Type> >(ptf);

    unitField_.rmap(tiptf.unitField_, addr);
}


template<class Type>
void analyticFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
	setCurField();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void analyticFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    unitField_.writeEntry("unitField", os);
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("waveVector")
        << waveVector_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
