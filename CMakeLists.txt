cmake_minimum_required(VERSION 3.23)
project(6)

set(CMAKE_CXX_STANDARD 14)

add_executable(6

        FEM/ErrorCompare/ErrorCompare.cpp FEM/ErrorCompare/ErrorCompare.h

        FEM/SolverMethod/SolveLinearEquations.cpp FEM/SolverMethod/SolveLinearEquations.h
        FEM/SolverMethod/IterativeMethod.cpp FEM/SolverMethod/IterativeMethod.h
        FEM/SolverMethod/DirectMethod.cpp FEM/SolverMethod/DirectMethod.h
        FEM/SolverMethod/SolveEigenValue.cpp FEM/SolverMethod/SolveEigenValue.h

        FEM/Tools/PrintTools.cpp FEM/Tools/PrintTools.h
        FEM/Tools/SymbolsOperat.h FEM/Tools/TypeConversion.h

        FEM/FEMCoreIdea/AssemblyMatrix.cpp FEM/FEMCoreIdea/AssemblyMatrix.h
        FEM/FEMCoreIdea/AssemblyVector.cpp FEM/FEMCoreIdea/AssemblyVector.h
        FEM/FEMCoreIdea/BoundaryInformation.cpp FEM/FEMCoreIdea/BoundaryInformation.h
        FEM/FEMCoreIdea/CoefficientFunction.cpp FEM/FEMCoreIdea/CoefficientFunction.h
        FEM/FEMCoreIdea/FeBasisFunction.cpp FEM/FEMCoreIdea/FeBasisFunction.h
        FEM/FEMCoreIdea/FEMSolver.cpp FEM/FEMCoreIdea/FEMSolver.h
        FEM/FEMCoreIdea/GenerateMatrixPbTb.cpp FEM/FEMCoreIdea/GenerateMatrixPbTb.h

        )
