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
 * \brief Test for the exercise on two-phase flow in fractured porous media.
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

// include the properties header
#include "properties.hh"

// Assuming the domain is saturated with CO2 before the injection!!!!
//#include "properties_pureco2.hh"


#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

// Define some types for this test  so that we can set them as properties below and
// reuse them again in the main function with only one single definition of them here
using MatrixTypeTag = Dumux::Properties::TTag::MatrixProblem;
using FractureTypeTag = Dumux::Properties::TTag::FractureProblem;
using MatrixFVGridGeometry = Dumux::GetPropType<MatrixTypeTag, Dumux::Properties::GridGeometry>;
using FractureFVGridGeometry = Dumux::GetPropType<FractureTypeTag, Dumux::Properties::GridGeometry>;
using TheMultiDomainTraits = Dumux::MultiDomainTraits<MatrixTypeTag, FractureTypeTag>;
using TheCouplingMapper = Dumux::FacetCouplingMapper<MatrixFVGridGeometry, FractureFVGridGeometry>;
using TheCouplingManager = Dumux::FacetCouplingManager<TheMultiDomainTraits, TheCouplingMapper>;

// set the coupling manager property in the sub-problems
namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::MatrixProblem> { using type = TheCouplingManager; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FractureProblem> { using type = TheCouplingManager; };

} // end namespace Properties
} // end namespace Dumux

// main program
int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // We use the grid manager from the facet coupling framework (see alias above main)
    // This requires the grids used to be passed as template arguments, where
    // they are assumed to be ordered in descending grid dimension. Thus,
    // we pass the matrix grid as first and the fracture grid as second argument.
    using MatrixGrid = GetPropType<MatrixTypeTag, Properties::Grid>;
    using FractureGrid = GetPropType<FractureTypeTag, Properties::Grid>;
    using GridManager = Dumux::FacetCouplingGridManager<MatrixGrid, FractureGrid>;
    GridManager gridManager;

    // try to create a grid (from the grid file in the input file)
    // Init() creates the grid from the grid file specified in the input file.
    // This works with a single grid file in which in addition to the matrix
    // (d-dimensional) elements also the (d-1)-dimensional elements are specified
    // which are interpreted as the fracture elements. See the .geo file in the
    // grids folder on how to create such grids using gmsh. Note that currently
    // only the Gmsh mesh format (.msh) is supported!
    gridManager.init();
    gridManager.loadBalance();

    // we compute on the leaf grid views (get them from grid manager)
    // the grid ids correspond to the order of the grids passed to the manager (see above)
    static constexpr std::size_t matrixGridId = 0;
    static constexpr std::size_t fractureGridId = 1;
    const auto& matrixGridView = gridManager.template grid<matrixGridId>().leafGridView();
    const auto& fractureGridView = gridManager.template grid<fractureGridId>().leafGridView();

    // create the finite volume grid geometries
    auto matrixFvGridGeometry = std::make_shared<MatrixFVGridGeometry>(matrixGridView);
    auto fractureFvGridGeometry = std::make_shared<FractureFVGridGeometry>(fractureGridView);
    matrixFvGridGeometry->update();
    fractureFvGridGeometry->update();

    // the problems (boundary/initial conditions etc)
    using MatrixProblem = GetPropType<MatrixTypeTag, Properties::Problem>;
    using FractureProblem = GetPropType<FractureTypeTag, Properties::Problem>;

    // pass the model parameter group to the spatial params so that they obtain the right
    // values from the input file since we use groups for matrix and fracture
    // We also want to use domain markers in this exercise. For this reason, we also pass
    // the grid data object from the grid manager to them, so that they have access to the
    // domain markers that were specified in the grid file.
    auto matrixGridData = gridManager.getGridData()->template getSubDomainGridData<matrixGridId>();
    auto matrixSpatialParams = std::make_shared<typename MatrixProblem::SpatialParams>(matrixFvGridGeometry, matrixGridData, "Matrix");
    auto matrixProblem = std::make_shared<MatrixProblem>(matrixFvGridGeometry, matrixSpatialParams, "Matrix");

    // extract domain height from the matrix problem and pass to fracture problem (needed for right initial pressure distribution)
    auto fractureGridData = gridManager.getGridData()->template getSubDomainGridData<fractureGridId>();
    auto fractureSpatialParams = std::make_shared<typename FractureProblem::SpatialParams>(fractureFvGridGeometry, fractureGridData, "Fracture");
    auto fractureProblem = std::make_shared<FractureProblem>(fractureFvGridGeometry, fractureSpatialParams, fractureGridData, "Fracture");

    // the solution vector
    using SolutionVector = typename TheMultiDomainTraits::SolutionVector;
    SolutionVector x, xOld;

    // The domain ids within the multi-domain framework.
    // They do not necessarily have to be the same as the grid ids
    // in case you have more subdomains involved. We domain ids
    // correspond to the order of the type tags passed to the multidomain
    // traits (see definition of the traits class at the beginning of this file)
    static const auto matrixDomainId = typename TheMultiDomainTraits::template SubDomain<0>::Index();
    static const auto fractureDomainId = typename TheMultiDomainTraits::template SubDomain<1>::Index();

    // resize the solution vector and write initial solution to it
    x[matrixDomainId].resize(matrixFvGridGeometry->numDofs());
    x[fractureDomainId].resize(fractureFvGridGeometry->numDofs());
    matrixProblem->applyInitialSolution(x[matrixDomainId]);
    fractureProblem->applyInitialSolution(x[fractureDomainId]);

    // instantiate the class holding the coupling maps between the domains
    // this needs the information on embeddings (connectivity between matrix
    // and fracture domain). This information is extracted directly from the
    // grid during file read and can therefore be obtained from the grid manager.
    const auto embeddings = gridManager.getEmbeddings();
    auto couplingMapper = std::make_shared<TheCouplingMapper>();
    couplingMapper->update(*matrixFvGridGeometry, *fractureFvGridGeometry, embeddings);

    // the coupling manager (needs the coupling mapper)
    auto couplingManager = std::make_shared<TheCouplingManager>();
    couplingManager->init(matrixProblem, fractureProblem, couplingMapper, x);

    // we have to set coupling manager pointer in sub-problems
    // they also have to be made accessible in them (see e.g. matrixproblem.hh)
    matrixProblem->setCouplingManager(couplingManager);
    fractureProblem->setCouplingManager(couplingManager);

    // the grid variables
    using MatrixGridVariables = GetPropType<MatrixTypeTag, Properties::GridVariables>;
    using FractureGridVariables = GetPropType<FractureTypeTag, Properties::GridVariables>;
    auto matrixGridVariables = std::make_shared<MatrixGridVariables>(matrixProblem, matrixFvGridGeometry);
    auto fractureGridVariables = std::make_shared<FractureGridVariables>(fractureProblem, fractureFvGridGeometry);
    xOld = x;
    matrixGridVariables->init(x[matrixDomainId]);
    fractureGridVariables->init(x[fractureDomainId]);

    // intialize the vtk output modules
    using MatrixVtkOutputModule = VtkOutputModule<MatrixGridVariables, GetPropType<MatrixTypeTag, Properties::SolutionVector>>;
    using FractureVtkOutputModule = VtkOutputModule<FractureGridVariables, GetPropType<FractureTypeTag, Properties::SolutionVector>>;
    MatrixVtkOutputModule matrixVtkWriter(*matrixGridVariables, x[matrixDomainId], matrixProblem->name());
    FractureVtkOutputModule fractureVtkWriter(*fractureGridVariables, x[fractureDomainId], fractureProblem->name());

    // Add model specific output fields
    using MatrixIOFields = GetPropType<MatrixTypeTag, Properties::IOFields>;
    using FractureIOFields = GetPropType<FractureTypeTag, Properties::IOFields>;
//    VtkOutputModule<MatrixGridVariables, SolutionVector> matrixVtkWriter(*matrixGridVariables, x[matrixDomainId], matrixProblem->name());
//    using VelocityOutput = GetPropType<MatrixTypeTag, Properties::VelocityOutput>;
//    matrixVtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*matrixGridVariables));
//	using VelocityOutput = GetPropType<FractureTypeTag, Properties::VelocityOutput>;
//	fractureVtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*fractureGridVariables));
    MatrixIOFields::initOutputModule(matrixVtkWriter);
    FractureIOFields::initOutputModule(fractureVtkWriter);

    // add domain markers to output
    std::vector<int> matrixDomainMarkers(matrixFvGridGeometry->gridView().size(0));
    for (const auto& element : elements(matrixFvGridGeometry->gridView()))
        matrixDomainMarkers[matrixFvGridGeometry->elementMapper().index(element)] = matrixProblem->spatialParams().getElementDomainMarker(element);
    matrixVtkWriter.addField(matrixDomainMarkers, "domainMarker");

    std::vector<int> fractureDomainMarkers(fractureFvGridGeometry->gridView().size(0));
    for (const auto& element : elements(fractureFvGridGeometry->gridView()))
        fractureDomainMarkers[fractureFvGridGeometry->elementMapper().index(element)] = fractureProblem->spatialParams().getElementDomainMarker(element);
    fractureVtkWriter.addField(fractureDomainMarkers, "domainMarker");

    // write out initial solution
    matrixProblem->updateVtkFields(x[matrixDomainId]);
    matrixProblem->addVtkFields(matrixVtkWriter);
    fractureProblem->updateVtkFields(x[fractureDomainId]);
    fractureProblem->addVtkFields(fractureVtkWriter); //!< Add problem specific output fields
    matrixVtkWriter.write(0.0);
    fractureVtkWriter.write(0.0);

//    using Velocity = Dune::FieldVector<double, 2>;
//    using VelocityVector = std::vector<Velocity>;

//	VelocityVector velocityMatrix(matrixFvGridGeometry->gridView().size(0));
//	VelocityVector velocityFracture(fractureFvGridGeometry->gridView().size(0));

//	matrixVtkWriter.addField(velocityMatrix, "velocity");
//	fractureVtkWriter.addField(velocityFracture, "velocity");

    // get some time loop parameters
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    const auto maxDt = getParam<double>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<double>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared< TimeLoop<double> >(/*startTime*/0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler for the coupled problem
    using Assembler = MultiDomainFVAssembler<TheMultiDomainTraits, TheCouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(matrixProblem, fractureProblem),
                                                  std::make_tuple(matrixFvGridGeometry, fractureFvGridGeometry),
                                                  std::make_tuple(matrixGridVariables, fractureGridVariables),
                                                  couplingManager,
                                                  timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, TheCouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        newtonSolver->solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        matrixGridVariables->advanceTimeStep();
        fractureGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

//        computeVelocities<MatrixTypeTag>(matrixDomainId, *assembler, *couplingManager, *matrixFvGridGeometry, *matrixGridVariables, x[matrixDomainId], velocityMatrix);
//		computeVelocities<FractureTypeTag>(fractureDomainId, *assembler, *couplingManager, *fractureFvGridGeometry, *fractureGridVariables, x[fractureDomainId], velocityFracture);

        matrixProblem->updateVtkFields(x[matrixDomainId]);
        fractureProblem->updateVtkFields(x[fractureDomainId]);

        // write vtk output
        matrixVtkWriter.write(timeLoop->time());
        fractureVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the Newton solver
        timeLoop->setTimeStepSize(newtonSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    // output some Newton statistics
    newtonSolver->report();

    // report time loop statistics
    timeLoop->finalize();

    // print dumux message to say goodbye
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}// end main
