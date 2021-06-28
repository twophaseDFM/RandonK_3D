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
 * \brief The properties file for exercise on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_PROPERTIES_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_PROPERTIES_HH

// Both sub-problems
// include the model we inherit from
#include <dumux/porousmediumflow/2p/model.hh>
// we want to simulate nitrogen gas transport in a water-saturated medium
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/brineco2.hh>

// Fracture sub-problem
// we use foam grid for the discretization of the fracture domain
// as this grid manager is able to represent network/surface grids
#include <dune/foamgrid/foamgrid.hh>
// we use a cell-centered finite volume scheme with tpfa here
#include <dumux/discretization/cctpfa.hh>
// the spatial parameters (permeabilities, material parameters etc.)
#include "fracturespatialparams.hh"
// the fracture sub-problem problem file
#include "fractureproblem.hh"

// Matrix sub-problem
// the spatial parameters (permeabilities, material parameters etc.)
#include "matrixspatialparams.hh"
// the matrix sub-problem problem file
#include "matrixproblem.hh"
#include "co2tables.hh"
// we use alu grid for the discretization of the matrix domain
#include <dune/alugrid/grid.hh>
// We are using the framework for models that consider coupling
// across the element facets of the bulk domain. This has some
// properties defined, which we have to inherit here. In this
// exercise we want to use a cell-centered finite volume scheme
// with tpfa.
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>

namespace Dumux::Properties {

// create the type tag node for the matrix and fracture sub-problems
namespace TTag {
struct MatrixProblem { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, TwoPNI>; };
struct FractureProblem { using InheritsFrom = std::tuple<TwoPNI, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type for matrix and fracture sub-domains
template<class TypeTag>
struct Grid<TypeTag, TTag::MatrixProblem> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::FractureProblem> { using type = Dune::FoamGrid<1, 2>; };

// Set the problem type for the matrix and fracture sub-domains
template<class TypeTag>
struct Problem<TypeTag, TTag::MatrixProblem> { using type = MatrixSubProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::FractureProblem> { using type = FractureSubProblem<TypeTag>; };

// set the spatial params for the matrix and fracture sub-domains
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MatrixProblem>
{
    using type = MatrixSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                      GetPropType<TypeTag, Properties::Scalar>>;
//									  GetPropType<TypeTag, Properties::GridVariables>>;
//	using type = MatrixSpatialParams< GetPropType<TypeTag>>;
};
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FractureProblem>
{
    using type = FractureSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
//										GetPropType<TypeTag, Properties::GridVariables>>;
//	using type = FractureSpatialParams<TypeTag>;
};

// the fluid system for the matrix and fracture sub-domains
template<class TypeTag>
//struct FluidSystem<TypeTag, TTag::MatrixProblem>
//{
//    using type = Dumux::FluidSystems::H2ON2< GetPropType<TypeTag, Properties::Scalar>,
//                                             FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true> >;
//};
struct FluidSystem<TypeTag, TTag::MatrixProblem>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        HeterogeneousCO2Tables::CO2Tables,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
};

template<class TypeTag>
//struct FluidSystem<TypeTag, TTag::FractureProblem>
//{
//    using type = Dumux::FluidSystems::H2ON2< GetPropType<TypeTag, Properties::Scalar>,
//                                             FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true> >;
//};
struct FluidSystem<TypeTag, TTag::FractureProblem>
{
    using type = FluidSystems::BrineCO2<GetPropType<TypeTag, Properties::Scalar>,
                                        HeterogeneousCO2Tables::CO2Tables,
                                        Components::TabulatedComponent<Components::H2O<GetPropType<TypeTag, Properties::Scalar>>>,
                                        FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
};

} // end namespace Dumux::Properties

#endif
