#ifndef BASICIO_H
#define BASICIO_H

/* file:                BasicIO.h
 * description:         functions for mesh loading and saving
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


#include <float.h>

namespace OFFIO{

template< class point_t > bool open( const std::string & filename ,
                                     std::vector< point_t > & vertices ,
                                     std::vector< std::vector< int > > & faces,
                                     bool convertToTriangles = true,
                                     bool randomize = false )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        return false;
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    vertices.resize(n_vertices);

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        vertices[v] = point_t( x , y , z );
    }


    for( int f = 0 ; f < n_faces ; ++f )
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;
        if( n_vertices_on_face == 3 )
        {
            int _v1 , _v2 , _v3;
            std::vector< int > _v;
            myfile >> _v1 >> _v2 >> _v3;
            _v.push_back( _v1 );
            _v.push_back( _v2 );
            _v.push_back( _v3 );
            faces.push_back( _v );
        }
        else if( n_vertices_on_face > 3 )
        {
            std::vector< int > vhandles;
            vhandles.resize(n_vertices_on_face);
            for( unsigned int i=0 ; i < n_vertices_on_face ; ++i )
                myfile >> vhandles[i];

            if( convertToTriangles )
            {
                unsigned int k=(randomize)?(rand()%vhandles.size()):0;
                for (unsigned int i=0;i<vhandles.size()-2;++i)
                {
                    std::vector< int > tri(3);
                    tri[0]=vhandles[(k+0)%vhandles.size()];
                    tri[1]=vhandles[(k+i+1)%vhandles.size()];
                    tri[2]=vhandles[(k+i+2)%vhandles.size()];
                    faces.push_back(tri);
                }
            }
            else
            {
                faces.push_back(vhandles);
            }
        }
        else
        {
            std::cout << "OFFIO::open error : Face number " << f << " has " << n_vertices_on_face << " vertices" << std::endl;
            myfile.close();
            return false;
        }
    }

    myfile.close();
    return true;
}


template< class point_t > bool save( const std::string & filename , std::vector< point_t > & vertices , std::vector< std::vector< int > > & faces )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = vertices.size() , n_faces = faces.size();
    myfile << n_vertices << " " << n_faces << " 0" << std::endl;

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( int f = 0 ; f < n_faces ; ++f )
    {
        myfile << faces[f].size();
        for( unsigned int vof = 0 ; vof < faces[f].size() ; ++vof )
            myfile << " " << faces[f][vof];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}


template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        return false;
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    vertices.resize(n_vertices);

    point_t bb( FLT_MAX , FLT_MAX , FLT_MAX );
    point_t BB( -FLT_MAX , -FLT_MAX , -FLT_MAX );

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        vertices[v] = point_t( x , y , z );
        bb = point_t::Min( bb , point_t( x , y , z ) );
        BB = point_t::Max( BB , point_t( x , y , z ) );
    }

    myfile.close();
}

}



#endif // BASICIO_H

