# First variables
BUILD_DIR = ./Build
CXX = g++
CXXFLAGS = -Wall -Wextra -Wpedantic -fPIC -O3 -pthread

.PHONY: clean

all: ${BUILD_DIR}/Combustion.so

${BUILD_DIR}/%.o: %.cpp %.hpp
	mkdir -p $(@D)
	${CXX} ${CXXFLAGS} -o $@ -c $<

${BUILD_DIR}/Math/Instantiation/Instantiation_double.o: ./Math/Instantiation/Instantiation_double.cpp\
		./Math/Vector.cpp\
		./Math/Vector.hpp\
		./Math/Matrix.cpp\
		./Math/Matrix.hpp\
		./Math/DualNumbers.cpp\
		./Math/DualNumbers.hpp
	mkdir -p ${BUILD_DIR}/Math/Instantiation
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/Math/Instantiation/Instantiation_double.o -c ./Math/Instantiation/Instantiation_double.cpp

${BUILD_DIR}/Math/Instantiation/Instantiation_complex.o: ./Math/Instantiation/Instantiation_complex.cpp\
		./Math/Vector.cpp\
		./Math/Vector.hpp
	mkdir -p ${BUILD_DIR}/Math/Instantiation
	${CXX} ${CXXFLAGS} -o ${BUILD_DIR}/Math/Instantiation/Instantiation_complex.o -c ./Math/Instantiation/Instantiation_complex.cpp

${BUILD_DIR}/Combustion.so: ${BUILD_DIR}/Math/Instantiation/Instantiation_double.o\
		${BUILD_DIR}/Math/Instantiation/Instantiation_complex.o\
		${BUILD_DIR}/Math/ODESolvers.o\
		${BUILD_DIR}/Math/Integration.o\
		${BUILD_DIR}/Math/AlgebraicSolvers.o\
		${BUILD_DIR}/Math/Interpolation.o\
		${BUILD_DIR}/Utilities/FileArray.o\
		${BUILD_DIR}/Utilities/FormatNumber.o\
		${BUILD_DIR}/Utilities/Progress.o\
		${BUILD_DIR}/CPGF/AffineSpace2d/Vector2d.o\
		${BUILD_DIR}/CPGF/AffineSpace2d/Point2d.o\
		${BUILD_DIR}/CPGF/PGFBasics/PGFConf.o\
		${BUILD_DIR}/CPGF/PGFBasics/Strokes2d.o\
		${BUILD_DIR}/CPGF/PGFBasics/Path2d.o\
		${BUILD_DIR}/CPGF/Text/Text.o\
		${BUILD_DIR}/CPGF/Objects2d/Object2d.o\
		${BUILD_DIR}/CPGF/Objects2d/BasicGeometries.o\
		${BUILD_DIR}/CPGF/Scene2d.o\
		${BUILD_DIR}/CPGF/Plot2d/Axis.o\
		${BUILD_DIR}/CPGF/Plot2d/DataPlot.o\
		${BUILD_DIR}/CPGF/Plot2d/Graphic.o\
		${BUILD_DIR}/CPGF/Plot2d/LinePlot.o\
		${BUILD_DIR}/Chemistry/Reaction.o\
		${BUILD_DIR}/Mesh/Cell.o\
		${BUILD_DIR}/Mesh/Mesh.o\
		${BUILD_DIR}/Solvers/Gas/ExactRiemannSolver.o\
		${BUILD_DIR}/Solvers/Gas/ExactSteadySolver.o\
		${BUILD_DIR}/Solvers/Gas/GasSolvers.o\
		${BUILD_DIR}/Solvers/Solid/SolidSolvers.o\
		${BUILD_DIR}/Solvers/Rocket/RocketSolver.o\
		${BUILD_DIR}/Simulation.o
	${CXX} -shared ${CXXFLAGS} -o ${BUILD_DIR}/Combustion.so\
		${BUILD_DIR}/Math/Instantiation/Instantiation_double.o\
		${BUILD_DIR}/Math/Instantiation/Instantiation_complex.o\
		${BUILD_DIR}/Math/ODESolvers.o\
		${BUILD_DIR}/Math/Integration.o\
		${BUILD_DIR}/Math/AlgebraicSolvers.o\
		${BUILD_DIR}/Math/Interpolation.o\
		${BUILD_DIR}/Utilities/FileArray.o\
		${BUILD_DIR}/Utilities/FormatNumber.o\
		${BUILD_DIR}/Utilities/Progress.o\
		${BUILD_DIR}/CPGF/AffineSpace2d/Vector2d.o\
		${BUILD_DIR}/CPGF/AffineSpace2d/Point2d.o\
		${BUILD_DIR}/CPGF/PGFBasics/PGFConf.o\
		${BUILD_DIR}/CPGF/PGFBasics/Strokes2d.o\
		${BUILD_DIR}/CPGF/PGFBasics/Path2d.o\
		${BUILD_DIR}/CPGF/Text/Text.o\
		${BUILD_DIR}/CPGF/Objects2d/Object2d.o\
		${BUILD_DIR}/CPGF/Objects2d/BasicGeometries.o\
		${BUILD_DIR}/CPGF/Scene2d.o\
		${BUILD_DIR}/CPGF/Plot2d/Axis.o\
		${BUILD_DIR}/CPGF/Plot2d/DataPlot.o\
		${BUILD_DIR}/CPGF/Plot2d/Graphic.o\
		${BUILD_DIR}/CPGF/Plot2d/LinePlot.o\
		${BUILD_DIR}/Chemistry/Reaction.o\
		${BUILD_DIR}/Mesh/Cell.o\
		${BUILD_DIR}/Mesh/Mesh.o\
		${BUILD_DIR}/Solvers/Gas/ExactRiemannSolver.o\
		${BUILD_DIR}/Solvers/Gas/ExactSteadySolver.o\
		${BUILD_DIR}/Solvers/Gas/GasSolvers.o\
		${BUILD_DIR}/Solvers/Solid/SolidSolvers.o\
		${BUILD_DIR}/Solvers/Rocket/RocketSolver.o\
		${BUILD_DIR}/Simulation.o

clean:
	rm -r ${BUILD_DIR}