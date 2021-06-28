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
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH

#include <dumux/io/grid/griddata.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include "fractureproblem.hh"
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial params the two-phase facet coupling test
 */
template< class FVGridGeometry, class Scalar >
//template<class TypeTag>
class FractureSpatialParams
: public FVSpatialParams< FVGridGeometry, Scalar, FractureSpatialParams<FVGridGeometry, Scalar> >
//: public FVSpatialParams <GetPropType<TypeTag, Properties::GridGeometry>,
//  	  	  	  	  	  	  GetPropType<TypeTag, Properties::Scalar>,
//						  FractureSpatialParams<TypeTag>>
{
//	using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
//	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
//	using ThisType = FractureSpatialParams<TypeTag>;
    using ThisType = FractureSpatialParams< FVGridGeometry, Scalar >;
    using ParentType = FVSpatialParams< FVGridGeometry, Scalar, ThisType >;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using IndexSet = typename GridView::IndexSet;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    //! export the type used for permeabilities
    using PermeabilityType = Scalar;

    //! the constructor
    FractureSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                          std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                          const std::string& paramGroup)
    : ParentType(fvGridGeometry)
    , gridDataPtr_(gridData)
    , pcKrSwCurve_("Fracture.SpatialParams")
    , barrierPcKrSwCurve_("Fracture.SpatialParams.Barrier")
    ,randomPermeability_(fvGridGeometry->gridView().size(dimWorld), 0.0)
    ,indexSet_(fvGridGeometry->gridView().indexSet())
    {
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        initRandomField(*fvGridGeometry);
    }

    //generate the geostatistic field using gstat
    void initRandomField(const FVGridGeometry& gg)
    {
        const auto& gridView = gg.gridView();
        const auto& elementMapper = gg.elementMapper();
        const auto gStatControlFile = getParam<std::string>("Gstat.ControlFile");
        const auto gStatInputFile = getParam<std::string>("Gstat.InputFile");
        const auto outputFilePrefix = getParam<std::string>("Gstat.OutputFilePrefix");

        // create random permeability object
        using RandomField = GstatRandomField<GridView, Scalar>;
        RandomField randomPermeabilityField(gridView, elementMapper);
        randomPermeabilityField.create(gStatControlFile,
                                       gStatInputFile,
                                       outputFilePrefix + ".dat",
                                       RandomField::FieldType::log10,
                                       true);
        randomPermeability_.resize(gridView.size(dimWorld), 0.0);

        // copy vector from the temporary gstat object
        randomPermeability_ = randomPermeabilityField.data();
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    template< class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
//    	return 1e-8;
    	return randomPermeability_[indexSet_.index(element)];
    }

    //! Return the porosity
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     *
     * \param element The current finite element
     * \param scv The sub-control volume
     * \param elemSol The current element solution
     */
//    template<class ElementSolution>
//    auto fluidMatrixInteraction(const Element& element,
//                                const SubControlVolume& scv,
//                                const ElementSolution& elemSol) const
//    {
//        return makeFluidMatrixInteraction(pcKrSwCurve_);
//    }

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

    Scalar temperature() const
    { return temperature_ ; }

    //! returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

    Scalar pressure(int phaseIdx) const
    {return pressure_ ;}

    const std::vector<Scalar>& getPermField() const
    { return randomPermeability_; }

//    const CouplingManager& couplingManager() const
//    { return *couplingManagerPtr_; }

private:
    //! pointer to the grid data (contains domain markers)
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
//    std::shared_ptr<CouplingManager> couplingManagerPtr_;

    Scalar porosity_;
    const PcKrSwCurve pcKrSwCurve_;
    const PcKrSwCurve barrierPcKrSwCurve_;
    Scalar pressure_;
    Scalar temperature_;
    std::vector<Scalar> randomPermeability_;
    const IndexSet& indexSet_;

};

} // end namespace Dumux

#endif
