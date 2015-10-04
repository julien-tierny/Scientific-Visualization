/* file:                main.cpp
 * description:         example for the computation of the derivatives of the Mean Value Coordinates in 3D.
 * author:              (C) Jean-Marc Thiery
 * date:                September 2013.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <MVCoordinates3D.h>
#include <point3.h>
#include <BasicIO.h>
int main( int argc , char** argv )
{
	// INPUT :
	std::string meshFileName;
	std::string cageFileName;

	int wrongUsage = 2;

	for( int i = 0 ; i < argc; ++i )
	{
		if( (std::string(argv[i]) == "-m")  &&  ( i+1 < argc ) )
		{
			meshFileName = std::string(argv[i+1]);
			--wrongUsage;
		}
		if( (std::string(argv[i]) == "-c")  &&  ( i+1 < argc ) )
		{
			cageFileName = std::string(argv[i+1]);
			--wrongUsage;
		}
	}

	if( wrongUsage )
	{
		std::cout << "Not enough arguments" << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << "./mvcTest -m <input surface OFF> -c <input cage OFF>" << std::endl;
		return 0;
	}

	// OUTPUT :
	const std::string & weightsFileName = "outputWeights.txt";
	const std::string & gradientsFileName = "outputGradients.txt";
	const std::string & HessiansFileName = "outputHessians.txt";

	// OPEN OUTPUT FILES :
	std::ofstream weightsFile;
    	weightsFile.open(weightsFileName.c_str());
    	if (!weightsFile.is_open())
    	{
        	std::cout << weightsFileName << " cannot be opened" << std::endl;
        	return false;
    	}

	std::ofstream gradientsFile;
    	gradientsFile.open(gradientsFileName.c_str());
    	if (!gradientsFile.is_open())
    	{
        	std::cout << gradientsFileName << " cannot be opened" << std::endl;
        	return false;
    	}

	std::ofstream HessiansFile;
    	HessiansFile.open(HessiansFileName.c_str());
    	if (!HessiansFile.is_open())
    	{
        	std::cout << HessiansFileName << " cannot be opened" << std::endl;
        	return false;
    	}

	// PARSE INPUT FILES AND CREATE STRUCTURE :
	std::vector< point3d > cageVertices , meshVertices , cageNormals;
	std::vector< std::vector< int > > cageTriangles;

	OFFIO::open( meshFileName , meshVertices );
	OFFIO::open( cageFileName , cageVertices , cageTriangles );

        cageNormals.resize( cageTriangles.size() );
        for( unsigned int t = 0 ; t < cageTriangles.size(); ++t )
        {
            cageNormals[ t ] = point3d::cross( cageVertices[ cageTriangles[t][1] ] - cageVertices[ cageTriangles[t][0] ] ,
                                                 cageVertices[ cageTriangles[t][2] ] - cageVertices[ cageTriangles[t][0] ] );
            cageNormals[ t ].normalize();
        }

	// COMPUTATION :
	unsigned int doItForThatManyVertices = 10;
	for( unsigned int v = 0 ; v < doItForThatManyVertices ; ++v )
	{
		std::vector< double > v_weights;
		std::vector< point3d > v_gradients;
		std::vector< mat33d > v_Hessians;
		MVCoordinates::MVC3D::computeCoordinatesAndDerivatives(
                    meshVertices[v],
                    cageTriangles,
                    cageVertices,
                    cageNormals,
                    v_weights,
                    v_gradients,
                    v_Hessians
                    );
		
		double wAccum = 0.0;
		point3d gAccum(0,0,0) , dPos(0,0,0);
		mat33d HAccum(0,0,0,0,0,0,0,0,0) , Jacobian(0,0,0,0,0,0,0,0,0) , Hx(0,0,0,0,0,0,0,0,0) , Hy(0,0,0,0,0,0,0,0,0) , Hz(0,0,0,0,0,0,0,0,0);
		for( unsigned int cv = 0 ; cv < cageVertices.size() ; ++cv )
		{
			wAccum += v_weights[cv];
			gAccum += v_gradients[cv];
			HAccum += v_Hessians[cv];
			dPos += v_weights[cv] * cageVertices[cv];
			Jacobian += mat33d::tensor(cageVertices[cv] , v_gradients[cv]);
			Hx += cageVertices[cv][0] * v_Hessians[cv];
			Hy += cageVertices[cv][1] * v_Hessians[cv];
			Hz += cageVertices[cv][2] * v_Hessians[cv];
			weightsFile << v_weights[cv] << " ";
			gradientsFile << v_gradients[cv] << " ";
			HessiansFile << v_Hessians[cv] << " ";
		}
		dPos -= meshVertices[v];
		std::cout << "POINT " << v << ":" << std::endl;
		std::cout << "sum of weights (should be 1):" << wAccum << std::endl;
		std::cout << "sum of gradients (should be 0):" << gAccum << std::endl;
		std::cout << "sum of Hessians (should be 0):" << std::endl << HAccum << std::endl;
		std::cout << "dPos (should be 0):" << dPos << std::endl;
		std::cout << "Jacobian (should be Id):" << std::endl << Jacobian << std::endl;
		std::cout << "Hx (should be 0):" << std::endl << Hx << std::endl;
		std::cout << "Hy (should be 0):" << std::endl << Hy << std::endl;
		std::cout << "Hz (should be 0):" << std::endl << Hz << std::endl << std::endl;
		weightsFile << std::endl;
		gradientsFile << std::endl;
		HessiansFile << std::endl;
	}

	weightsFile.close();
	gradientsFile.close();
	HessiansFile.close();

    	return 1;
}

