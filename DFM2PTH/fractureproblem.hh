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
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
  */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/io/grid/griddata.hh>
#include <dune/common/indices.hh>
#include "fracturespatialparams.hh"
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class FractureSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    static constexpr int dimWorld = GridView::dimensionworld;


    // some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,

        //! Equation indices
        contiCO2EqIdx = Indices::conti0EqIdx + FluidSystem::CO2Idx,
		contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::BrineIdx,
        energyEqIdx = Indices::energyEqIdx,

        wPhaseIdx = FluidSystem::BrineIdx,
        nPhaseIdx = FluidSystem::CO2Idx,
    };

public:
    //! The constructor
    FractureSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                       std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
					   std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                       const std::string& paramGroup)
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , gridDataPtr_(gridData)
    , aperture1_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture1"))
    , aperture2_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture2"))
    , aperture3_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture3"))
    , aperture4_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture4"))
    , aperture5_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture5"))
    , IsPureCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsPureCO2"))
    , kn_(getParamFromGroup<Scalar>(paramGroup, "Problem.Stiffness"))
    , wte_(getParamFromGroup<Scalar>(paramGroup, "Problem.WettingPhaseThermalExpansionCoefficient"))
    , nte_(getParamFromGroup<Scalar>(paramGroup, "Problem.NonWettingPhaseThermalExpansionCoefficient"))
    , a_(getParamFromGroup<Scalar>(paramGroup, "Problem.a"))
    , b_(getParamFromGroup<Scalar>(paramGroup, "Problem.b"))
//    , InjectionRate_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionRate"))
//    , InjectionTemperature_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionTemperature"))
//    , IsInjectCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsInjectCO2"))
    {
        // initialize the fluid system, i.e. the tabulation
        // of water properties. Use the default p/T ranges.
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        FluidSystem::init();

        using PermeabilityType = Scalar;
    }

    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        vtk.addField(deltaT_, "deltaT");
        vtk.addField(overPressure_, "deltaP");

    }

    void updateVtkFields(const SolutionVector& curSol)
    {
    	int dofCodim = 0;
    	const auto& gridView = this->gridGeometry().gridView();
    	deltaT_.resize(gridView.size(dofCodim));
    	overPressure_.resize(gridView.size(dofCodim));

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
//
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                const GlobalPosition& globalPos = scv.center();
//                VolumeVariables volVars;
                const auto& priVars = elemSol[scv.localDofIndex()];
                const auto initialValues = initialAtPos(globalPos);
//                volVars.update(elemSol, *this, element, scv);
                deltaT_[dofIdxGlobal] = priVars[temperatureIdx] - initialValues[temperatureIdx];
                overPressure_[dofIdxGlobal] = priVars[pressureIdx] - initialValues[pressureIdx];
            }
        }
    }

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
            values.setAllDirichlet();

        return values;
    }

    //! Evaluate the source term in a sub-control volume of an element
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain using the function in the coupling manager
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);

        // these sources are in kg/s, divide by volume and extrusion to have it in kg/s/m³
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    template< class ElementSolution >
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
//	    FluidState fs;
//        VolumeVariables volVars;
//		const auto peff_ = volVars.saturation(nPhaseIdx) * volVars.pressure(nPhaseIdx) + volVars.saturation(wPhaseIdx)* volVars.pressure(wPhaseIdx);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto fmi = this->spatialParams().fluidMatrixInteraction(element, scv, elemSol);
        const auto peff_ = priVars[saturationIdx] * (priVars[pressureIdx] + fmi.pc(1 - priVars[saturationIdx])) + (1 - priVars[saturationIdx]) * priVars[pressureIdx];

		const GlobalPosition& globalPos = scv.center();
        const auto initialValues = initialAtPos(globalPos);
        const auto ThermalExpan = (priVars[temperatureIdx] - initialValues[temperatureIdx]) * (priVars[saturationIdx] * nte_ + (1 - priVars[saturationIdx]) * wte_);
        const auto deltaP = peff_ - initialValues[pressureIdx];
//        const auto deltaT_ = volVars.temperature() - initialValues[temperatureIdx];

		if (getElementDomainMarker(element) == 1)
			return aperture1_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture1_ ;
		else if (getElementDomainMarker(element) == 2)
			return aperture2_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture2_;
		else if (getElementDomainMarker(element) == 3)
			return aperture3_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture3_;
		else if (getElementDomainMarker(element) == 4)
			return aperture4_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture4_;
		else
        	return aperture5_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture5_;
//    	return 1e-3;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // For the grid used here, the height of the domain is equal
        // to the maximum y-coordinate
        const auto domainHeight = this->gridGeometry().bBoxMax()[dimWorld-1] + 6000;
        // we assume a constant water density of 1000 for initial conditions!
        const auto& g = this->spatialParams().gravity(globalPos);
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 - (domainHeight - globalPos[dimWorld-1])*densityW*g[dimWorld-1];
        values[temperatureIdx] = 283.0 + (domainHeight - globalPos[dimWorld-1])*0.03;
        if (IsPureCO2_)
        	{values[saturationIdx] = 1.0 - eps_;}
        else
        	{values[saturationIdx] = 0.0;}

        return values;
    }

//    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
//    {
//        NumEqVector values(0.0);
//
//        if (globalPos[dimWorld-1] < 75 + eps_ && globalPos[dimWorld-1] > 25 - eps_ && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
//        {
//            // compute enthalpy flux associated with this injection [(J/(kg*s)]
//            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
//            FluidState fs;
//
//            const auto initialValues = initialAtPos(globalPos);
//            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
//            fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
//            fs.setTemperature(wPhaseIdx, InjectionTemperature_);
//            fs.setTemperature(nPhaseIdx, InjectionTemperature_);
//
//            // energy flux is mass flux times specific enthalpy
//            if (IsInjectCO2_)
//                // inject air. negative values mean injection
//			{
//            	values[contiCO2EqIdx] = -InjectionRate_ * 700; // kg/(s*m^2) flow rate * density
//            	values[energyEqIdx] = values[contiCO2EqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
//			}
//			else
//			{
//                values[contiH2OEqIdx] = -InjectionRate_ * 1000; // kg/(s*m^2) flow rate * density
//            	values[energyEqIdx] = values[contiH2OEqIdx]*FluidSystem::enthalpy(fs, wPhaseIdx);
//			}
//        }
//
//        return values;
//    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return 283.15; /*10°*/ }

//    Scalar pressure(int PhaseIdx) const
//    {return pressure_ ;}

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
    Scalar aperture1_,aperture2_,aperture3_,aperture4_,aperture5_;
    bool IsPureCO2_;
//    bool IsInjectCO2_;
    static constexpr Scalar eps_ = 1e-7;
    Scalar kn_, nte_, wte_;
    Scalar temperature_;
    Scalar a_, b_;
    std::vector<Scalar> deltaT_ , overPressure_;
//    Scalar InjectionTemperature_;
//    Scalar InjectionRate_;
};

} // end namespace Dumux

#endif
