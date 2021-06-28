// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial parameters for the fracture sub-domain in the exercise
 *        on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_SPATIALPARAMS_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_SPATIALPARAMS_HH

#include <dumux/io/grid/griddata.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial params the two-phase facet coupling test
 */
template< class FVGridGeometry, class Scalar >
class MatrixSpatialParams
: public FVSpatialParams< FVGridGeometry, Scalar, MatrixSpatialParams<FVGridGeometry, Scalar> >
{
    using ThisType = MatrixSpatialParams< FVGridGeometry, Scalar >;
    using ParentType = FVSpatialParams< FVGridGeometry, Scalar, ThisType >;

    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    //! export the type used for permeabilities
    using PermeabilityType = Scalar;

    //! the constructor
    MatrixSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                        const std::string& paramGroup)
    : ParentType(fvGridGeometry)
    , gridDataPtr_(gridData)
    , pcKrSwCurve_("Matrix.SpatialParams")
    {
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        permeabilityRock_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.PermeabilityRock");

    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
//    	if (globalPos[0] < 0 + eps_)
//    		return permeabilityWell_ ;
//    	else
		return permeabilityRock_;
    }

    //! Return the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {   
        return makeFluidMatrixInteraction(pcKrSwCurve_); 
    }

    //! Water is the wetting phase
    template< class FluidSystem >
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        // we set water as the wetting phase here
        // which is phase0Idx in the H2oN2 fluid system
        return FluidSystem::phase0Idx;
    }

    //! returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    //! pointer to the grid data (contains domain markers)
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;

    Scalar porosity_;
    PermeabilityType permeabilityRock_;
    const PcKrSwCurve pcKrSwCurve_;
    static constexpr Scalar eps_ = 1e-7 ;
};

} // end namespace Dumux

#endif
