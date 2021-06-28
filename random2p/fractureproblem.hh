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
#include <cmath>
#include <algorithm>

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
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;


    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

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
    , IsPureCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsPureCO2"))
    , kn_(getParamFromGroup<Scalar>(paramGroup, "Problem.Stiffness"))
    , wte_(getParamFromGroup<Scalar>(paramGroup, "Problem.WettingPhaseThermalExpansionCoefficient"))
    , nte_(getParamFromGroup<Scalar>(paramGroup, "Problem.NonWettingPhaseThermalExpansionCoefficient"))
    , a_(getParamFromGroup<Scalar>(paramGroup, "Problem.a"))
    , b_(getParamFromGroup<Scalar>(paramGroup, "Problem.b"))
    , InjectionRate_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionRate"))
    , InjectionTemperature_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionTemperature"))
    , IsInjectCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsInjectCO2"))
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
        vtk.addField(this->spatialParams().getPermField(), "K");

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
        values.setAllDirichlet();

        if (isInjectionWell(globalPos) || isProductionWell(globalPos))
        	values.setAllNeumann();

//        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
//            values.setAllDirichlet();

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
        return source;//            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
        //            FluidState fs;
        //            const Scalar outletT_ = min(volVars.temperature(),initialValues[temperatureIdx]);
        //            fs.setPressure(wPhaseIdx, volVars.pressure(wPhaseIdx));
        //            fs.setPressure(nPhaseIdx, volVars.pressure(nPhaseIdx)); // assume pressure equality here
        //            fs.setTemperature(wPhaseIdx, outletT_);
        //            fs.setTemperature(nPhaseIdx, outletT_);
    }

    template< class ElementSolution >
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
    	return 1e-3;
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

    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& globalPos = scvf.ipGlobal();
        const auto& fluidMatrixInteraction = this->spatialParams().fluidMatrixInteractionAtPos(globalPos);
        const auto initialValues = initialAtPos(globalPos);

        if (isInjectionWell(globalPos))
        {
            // compute enthalpy flux associated with this injection [(J/(kg*s)]
            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
            FluidState fs;

            const auto pn = initialValues[pressureIdx] + fluidMatrixInteraction.endPointPc();
            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
            fs.setPressure(nPhaseIdx, pn); // assume pressure equality here
            fs.setTemperature(wPhaseIdx, InjectionTemperature_);
            fs.setTemperature(nPhaseIdx, InjectionTemperature_);

            // energy flux is mass flux times specific enthalpy
            if (IsInjectCO2_)
                // inject air. negative values mean injection
			{
            	values[contiCO2EqIdx] = -InjectionRate_ * 700; // kg/(s*m^2) flow rate * density
            	values[energyEqIdx] = values[contiCO2EqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
			}
			else
			{
                values[contiH2OEqIdx] = -InjectionRate_ * 1000; // kg/(s*m^2) flow rate * density
            	values[energyEqIdx] = values[contiH2OEqIdx]*FluidSystem::enthalpy(fs, wPhaseIdx);
			}
        }

        if (isProductionWell(globalPos))
        {
            // compute enthalpy flux associated with this injection [(J/(kg*s)]
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            using std::max;
            using std::min;
            const Scalar sneff_ = max(volVars.saturation(nPhaseIdx),0.0);
//            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
//            FluidState fs;
//            const Scalar outletT_ = min(volVars.temperature(),initialValues[temperatureIdx]);
//            fs.setPressure(wPhaseIdx, volVars.pressure(wPhaseIdx));
//            fs.setPressure(nPhaseIdx, volVars.pressure(nPhaseIdx)); // assume pressure equality here
//            fs.setTemperature(wPhaseIdx, outletT_);
//            fs.setTemperature(nPhaseIdx, outletT_);

            // energy flux is mass flux times specific enthalpy
			values[contiCO2EqIdx] = sneff_ * InjectionRate_ * 700; // kg/(s*m^2) flow rate * density
			values[contiH2OEqIdx] = volVars.saturation(wPhaseIdx) * InjectionRate_ * 1000;
			values[energyEqIdx] = values[contiCO2EqIdx] * FluidSystem::enthalpy(volVars.fluidState(),nPhaseIdx) +
					              values[contiH2OEqIdx] * FluidSystem::enthalpy(volVars.fluidState(),wPhaseIdx);
        }


        return values;
    }

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
    bool IsInjectCO2_;
    static constexpr Scalar eps_ = 1e-7;
    Scalar kn_, nte_, wte_;
    Scalar temperature_;
    Scalar a_, b_;
    std::vector<Scalar> deltaT_ , overPressure_;
    Scalar InjectionTemperature_;
    Scalar InjectionRate_;

    bool isInjectionWell (const GlobalPosition& globalPos) const
    {
    	return globalPos[1] < 55 + eps_ && globalPos[1] > 45 - eps_ && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool isProductionWell (const GlobalPosition& globalPos) const
    {
    	return globalPos[1] < 55 + eps_ && globalPos[1] > 45 - eps_ && globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }
};

} // end namespace Dumux

#endif
