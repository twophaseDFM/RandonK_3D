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
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH

// we need this in this test in order to define the domain
// id of the fracture problem (see function interiorBoundaryTypes())
#include <dune/common/indices.hh>

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class MatrixSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = typename GridVariables::Scalar;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

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

        //! Phase indices
        wPhaseIdx = FluidSystem::BrineIdx,
        nPhaseIdx = FluidSystem::CO2Idx,
    };

public:
    //! The constructor
    MatrixSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                     const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
//    , boundaryOverPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryOverPressure"))
//    , boundarySaturation_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundarySaturation"))
    , InjectionRate_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionRate"))
    , InjectionTemperature_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionTemperature"))
    , IsInjectCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsInjectCO2"))
    , IsPureCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsPureCO2"))
    {
        // initialize the fluid system, i.e. the tabulation
        // of water properties. Use the default p/T ranges.
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        FluidSystem::init();
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
//
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

        // we consider buoancy-driven upwards migration of nitrogen and set
        // Dirichlet BCs on the top and bottom boundary
        values.setAllNeumann();
//        if (globalPos[1] < 1e-6 || globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
		if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_ && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_ )
            values.setAllDirichlet();

        return values;
    }

    //! Specifies the type of interior boundary condition to be used on a sub-control volume face
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // initialize values with the initial conditions
        auto values = initialAtPos(globalPos);

        return values;
    }

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

    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        return values;
    }

//    void addPointSources(std::vector<PointSource>& pointSources) const
//    {
//        // inject 10 kg/s water at origin;
//        pointSources.push_back(PointSource(GlobalPosition(25,50,50), {10}));
//    }
//
//    NumEqVector source(const Element &element,
//                       const FVElementGeometry& fvGeometry,
//                       const ElementVolumeVariables& elemVolVars,
//                       const SubControlVolume &scv) const
//    {
//        NumEqVector source(0.0);
//
//		using FluidState = GetPropType<TypeTag, Properties::FluidState>;
//		FluidState fs;
//
//        const auto& globalPos = scv.center();
//        const auto& volVars = elemVolVars[scv];
//
//		const auto initialValues = initialAtPos(globalPos);
//		fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
//		fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
//		fs.setTemperature(wPhaseIdx, InjectionTemperature_);
//		fs.setTemperature(nPhaseIdx, InjectionTemperature_);
//
//        if (isInjectionWell(globalPos))
//            source[contiCO2EqIdx] = -InjectionRate_;
//        	source[energyEqIdx] = InjectionRate_ * FluidSystem::enthalpy(fs,nPhaseIdx);
//
//        return source;
//    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return temperature_; }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;

    Scalar boundaryOverPressure_;
    Scalar boundarySaturation_;
    Scalar InjectionTemperature_;
    Scalar InjectionRate_;
    bool IsInjectCO2_;
    bool IsPureCO2_;
    Scalar temperature_;
    Scalar saturation_;
    static constexpr Scalar eps_ = 1e-7;
    std::vector<Scalar> overPressure_, deltaT_;

    bool isInjectionWell (const GlobalPosition& globalPos) const
    {
    	return globalPos[0] < 30 && globalPos[0] > 20 && globalPos[1] < 51 && globalPos[1] > 49 && globalPos[2] < 51 && globalPos[2] > 49;
    }

    bool isProductionWell (const GlobalPosition& globalPos) const
    {
    	return globalPos[0] < 80 && globalPos[0] > 70 && globalPos[1] < 51 && globalPos[1] > 49 && globalPos[2] < 51 && globalPos[2] > 49;
    }
};

} // end namespace Dumux

#endif
