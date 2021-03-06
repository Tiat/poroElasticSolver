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

Class
    Foam::analyticFixedValueFvPatchField

Description
    Foam::analyticFixedValueFvPatchField is a modified version of 
    Foam::oscillatingFixedValueFvPatchField, which provides the option to set 
    the value of a field on a patch to some analytical expression which may 
    both depend on spatial coordinates and time. 

SourceFiles
    analyticFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef analyticFixedValueFvPatchField_H
#define analyticFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class analyticFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class analyticFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> unitField_;

        //- Amplitude
        scalar amplitude_;

        //- Frequency
        scalar frequency_;

        //- Frequency
        vector waveVector_;

        //- Current time index
        label curTimeIndex_;


    // Private member functions

        //- Set current field values on patch
        void setCurField();

public:

    //- Runtime type information
    TypeName("analyticFixedValue");


    // Constructors

        //- Construct from patch and internal field
        analyticFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        analyticFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given analyticFixedValueFvPatchField
        //  onto a new patch
        analyticFixedValueFvPatchField
        (
            const analyticFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        analyticFixedValueFvPatchField
        (
            const analyticFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new analyticFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        analyticFixedValueFvPatchField
        (
            const analyticFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new analyticFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& unitField() const
            {
                return unitField_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& unitField()
            {
                return unitField_;
            }

            //- Return amplitude
            scalar amplitude() const
            {
                return amplitude_;
            }

            scalar& amplitude()
            {
                return amplitude_;
            }

            //- Return frequency
            scalar frequency() const
            {
                return frequency_;
            }

            scalar& frequency()
            {
                return frequency_;
            }

            //- Return wave vector
            scalar waveVector() const
            {
                return waveVector_;
            }

            scalar& waveVector()
            {
                return waveVector_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "analyticFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
