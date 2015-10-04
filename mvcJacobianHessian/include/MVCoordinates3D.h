#ifndef MVCOORDINATES3D_H
#define MVCOORDINATES3D_H

/* file:                MVCoordinates3D.h
 * description:         functions computing the derivatives of the Mean Value Coordinates in 3D.
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



#include "MVC_Equiv.h"

#include <vector>
#include <cmath>

#include <cassert>



namespace MVCoordinates{
namespace MVC3D{

// MVC :
// Code from "Mean Value Coordinates for Closed Triangular Meshes" Schaeffer Siggraph 2005
template< class point_t >
void computeCoordinatesOriginalCode(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights)
{
    typedef typename point_t::type_t    T;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    T epsilon = 0.00000001;

    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    std::vector< T > w_weights(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    T w[3];
    T c[3];
    T s[3];


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            vid[i] =  cage_triangles[t][i];
     //       std::cout << "vid " << vid[i] << std::endl;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();
      //      std::cout << "l " << l[i] << std::endl;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
     //       std::cout << "theta " << theta[i] << std::endl;
        }

        T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;

        if( M_PI - h < epsilon )
        {
            // eta is on the triangle t , use 2d barycentric coordinates :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
                w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];

            sumWeights = w[0] + w[1] + w[2];

            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[ vid[0] ] = w[0] / sumWeights;
            weights[ vid[1] ] = w[1] / sumWeights;
            weights[ vid[2] ] = w[2] / sumWeights;
            return;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            c[ i ] = ( 2.0 * sin(h) * sin(h - theta[ i ]) ) / ( sin(theta[ (i+1) % 3 ]) * sin(theta[ (i+2) % 3 ]) ) - 1.0;
   //         std::cout << "c " << c[i] << std::endl;
        }

        T sign_Basis_u0u1u2 = 1;
        if( point_t::dot( point_t::cross(u[vid[0]] , u[vid[1]]) , u[vid[2]] ) < 0.0 )
            sign_Basis_u0u1u2 = -1;

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            s[ i ] = sign_Basis_u0u1u2 * sqrt( std::max( (T)(0.0) , (T)(1.0) - c[ i ] * c[ i ] ) );
    //        std::cout << "s " << s[i] << std::endl;
        }

        if( std::abs( s[0] ) < epsilon   ||   std::abs( s[1] ) < epsilon   ||   std::abs( s[2] ) < epsilon )
        {
            // eta is on the same plane, outside t  ->  ignore triangle t :
            continue;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            w[ i ] = ( theta[ i ] - c[ (i+1)% 3 ]*theta[ (i+2) % 3 ] - c[ (i+2) % 3 ]*theta[ (i+1) % 3 ] ) / ( 2.0 * d[ vid[i] ] * sin( theta[ (i+1) % 3 ] ) * s[ (i+2) % 3 ] );

     //      std::cout << "w " << w[i] << std::endl;
        }

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        weights[v]  = w_weights[v] / sumWeights;
    }
}









template< class point_t >
void computeCoordinatesOriginalCode(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< typename point_t::type_t > & w_weights)
{
    typedef typename point_t::type_t    T;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    T epsilon = 0.00000001;

    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    w_weights.clear();
    w_weights.resize(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    T w[3];
    T c[3];
    T s[3];


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            vid[i] =  cage_triangles[t][i];

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            theta[i] = 2.0 * asin( l[i] / 2.0 );

        T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;

        if( M_PI - h < epsilon )
        {
            // eta is on the triangle t , use 2d barycentric coordinates :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
                w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];

            sumWeights = w[0] + w[1] + w[2];

            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[ vid[0] ] = w[0] / sumWeights;
            weights[ vid[1] ] = w[1] / sumWeights;
            weights[ vid[2] ] = w[2] / sumWeights;
            return;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            c[ i ] = ( 2.0 * sin(h) * sin(h - theta[ i ]) ) / ( sin(theta[ (i+1) % 3 ]) * sin(theta[ (i+2) % 3 ]) ) - 1.0;

        T sign_Basis_u0u1u2 = 1;
        if( point_t::dot( point_t::cross(u[vid[0]] , u[vid[1]]) , u[vid[2]] ) < 0.0 )
            sign_Basis_u0u1u2 = -1;

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            s[ i ] = sign_Basis_u0u1u2 * sqrt( max( 0.0 , 1.0 - c[ i ] * c[ i ] ) );

        if( abs( s[0] ) < epsilon   ||   abs( s[1] ) < epsilon   ||   abs( s[2] ) < epsilon )
        {
            // eta is on the same plane, outside t  ->  ignore triangle t :
            continue;
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            w[ i ] = ( theta[ i ] - c[ (i+1)% 3 ]*theta[ (i+2) % 3 ] - c[ (i+2) % 3 ]*theta[ (i+1) % 3 ] ) / ( 2.0 * d[ vid[i] ] * sin( theta[ (i+1) % 3 ] ) * s[ (i+2) % 3 ] );

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        weights[v]  = w_weights[v] / sumWeights;
    }
}









template< class point_t >
void computeCoordinatesSimpleCode(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights)
{
    typedef typename point_t::type_t    T;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    T epsilon = 0.000000001;

    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    std::vector< T > w_weights(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    T w[3];


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            vid[i] =  cage_triangles[t][i];

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
        }
        // test in original MVC paper: (they test if one angle psi is close to 0: it is "distance sensitive" in the sense that it does not
        // relate directly to the distance to the support plane of the triangle, and the farther away you go from the triangle, the worse it is)

        // simple test we suggest:
        // the determinant of the basis is 2*area(T)*d( eta , support(T) ), we can directly test for the distance to support plane of the triangle to be minimum
        T determinant = point_t::dot( cage_vertices[vid[0]] - eta , point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
        T sqrdist = determinant*determinant / (4 * point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ).sqrnorm() );
        T dist = sqrt( (T)sqrdist );

        if( dist < epsilon  &&  false )
        {
            // eta is on the triangle t , use 2d barycentric coordinates :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
            {
                w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
            }
            sumWeights = w[0] + w[1] + w[2];

            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[ vid[0] ] = w[0] / sumWeights;
            weights[ vid[1] ] = w[1] / sumWeights;
            weights[ vid[2] ] = w[2] / sumWeights;
            return;
        }

        point_t pt[3] , N[3];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            w[i] = 0.0;
            for( unsigned int j = 0 ; j <= 2 ; ++j )
            {
                w[i] += theta[j] * point_t::dot( N[i] , N[j] ) / ( 2 * N[j].norm() );
            }
            w[i] /= determinant;
        }

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;
}



template< class point_t >
void computeCoordinatesSimpleCode(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< typename point_t::type_t > & w_weights)
{
    typedef typename point_t::type_t    T;

    T epsilon = 0.000000001;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    w_weights.clear( );
    w_weights.resize(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    T w[3];


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            vid[i] =  cage_triangles[t][i];

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
        }

        T determinant = point_t::dot( cage_vertices[vid[0]] - eta , point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
        point_t pt[3] , N[3];

        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            w[i] = 0.0;
            for( unsigned int j = 0 ; j <= 2 ; ++j )
            {
                w[i] += theta[j] * point_t::dot( N[i] , N[j] ) / ( 2.0 * N[j].norm() );
            }
            w[i] /= determinant;
        }

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;
}








template< class point_t >
inline
void computeCoordinates(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights)
{
    computeCoordinatesOriginalCode(eta,cage_triangles,cage_vertices,cage_normals,weights);
//    computeCoordinatesSimpleCode(eta,cage_triangles,cage_vertices,cage_normals,weights);
}


template< class point_t >
inline
void computeCoordinates(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< typename point_t::type_t > & w_weights)// unnormalized weights
{
//    computeCoordinatesOriginalCode(eta,cage_triangles,cage_vertices,cage_normals,weights,w_weights);
    computeCoordinatesSimpleCode(eta,cage_triangles,cage_vertices,cage_normals,weights,w_weights);
}

















template< class point_t , class mat_t >
void computeCoordinatesAndDerivatives(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< point_t > & gradients ,
    std::vector< mat_t > & Hessians)
{
    typedef typename point_t::type_t    T;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    T epsilon = 0.000000001;

    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    std::vector< T > w_weights(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    std::vector< T > t_theta( 3*n_triangles , 0.0 );
    T w[3];
    T c[3];
    T s[3];
    std::vector< bool > point_is_on_triangle_support(n_triangles, false);
    std::vector< T > sum_W_on_triangles(n_triangles , 0.0);


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            vid[i] =  cage_triangles[t][i];

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
            // maybe change the evaluation of theta, maybe keep sin_theta and express the Taylor expansion formula in terms of sin_theta
            t_theta[ 3*t + i ] = theta[i];
        }

        T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;

        if( M_PI - h < epsilon )
        {
            // eta is on the triangle t , use 2d barycentric coordinates :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
            {
                w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
            }
            sumWeights = w[0] + w[1] + w[2];

            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[ vid[0] ] = w[0] / sumWeights;
            weights[ vid[1] ] = w[1] / sumWeights;
            weights[ vid[2] ] = w[2] / sumWeights;
            return;
        }

        // cos( psi_i ): corresponds to the dot product of the normals of the two adjacent spherical triangles.
        // they say in the paper that it is more robust to compute it that way...
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            c[ i ] = ( 2.0 * sin(h) * sin(h - theta[ i ]) ) / ( sin(theta[ (i+1) % 3 ]) * sin(theta[ (i+2) % 3 ]) ) - 1.0;

        T sign_Basis_u0u1u2 = 1;
        if( point_t::dot( point_t::cross(u[vid[0]] , u[vid[1]]) , u[vid[2]]) < 0.0 )
            sign_Basis_u0u1u2 = -1;

        // sin( psi_i ):
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            s[ i ] = sign_Basis_u0u1u2 * sqrt( std::max( (T)(0.0) , (T)(1.0) - c[ i ] * c[ i ] ) );

        // test in original MVC paper: (they test if one angle psi is close to 0: it is "distance sensitive" in the sense that it does not
        // relate directly to the distance to the support plane of the triangle, and the farther away you go from the triangle, the worse it is)

        // simple test we suggest:
        // the determinant of the basis is 2*area(T)*d( eta , support(T) ), we can directly test for the distance to support plane of the triangle to be minimum
        T determinant = point_t::dot( cage_vertices[vid[0]] - eta , point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
        T sqrdist = determinant*determinant / (4 * point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ).sqrnorm() );
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
        {
            point_is_on_triangle_support[t] = true;
            sum_W_on_triangles[t] = 0.0;
            continue;
        }
        // Note that epsilon should be multiplied by the diagonal of the cage model, in order to be correct
        // we cannot simply transform the model in the unit box: it does not change weights, but it scales Jacobians (linearly) and Hessians (quadratically)

        // maybe this should be reconsidered as well, as it does not increase at all the precision in the computations.
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            w[ i ] = ( theta[ i ] - c[ (i+1)% 3 ]*theta[ (i+2) % 3 ] - c[ (i+2) % 3 ]*theta[ (i+1) % 3 ] ) / ( 2.0 * d[ vid[i] ] * sin( theta[ (i+1) % 3 ] ) * s[ (i+2) % 3 ] );

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];

        sum_W_on_triangles[t] = ( w[0] + w[1] + w[2] );
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;

    // Computation of weights is done.



    std::vector< point_t > w_gradients( n_vertices , point_t(0,0,0) );
    std::vector< point_t > sum_gradW_on_triangles( n_triangles , point_t(0,0,0) );
    std::vector< mat_t > w_Hessians( n_vertices , mat_t(0,0,0,0,0,0,0,0,0) );

    // Code for Derivatives :
    point_t pt[3];  point_t N[3];
    T _eq1[3];      T _eq2[3];      T _eq3[3];      T _eq4[3];      T _eq5[3];      T _eq6[3];      T _eq7[3];      T _eq8[3];      T _eq9[3];      T _cos[3];
    T dt[3];        T dt_2[3];      T dt_3[3];      T dt_4[3];      T dt_5[3];
    T uu[3][3];
    T NN[3][3];
    point_t JNN[3][3];
    mat_t JN[3];






    const point_t  & NT = point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] );
    T NT_sqrnorm = point_t::dot( NT,NT );
    // ! NT_sqrnorm = 4 |T|^2 = NT_norm * 2 |T|
    const point_t & _nT_NT = - NT / (NT_sqrnorm);
    // _nT_NT = nT / (-2|T|)

    for( unsigned int i = 0 ; i < 3 ; ++i )
    {
        _eq1[i] = MVCD__EQUIV::eq1( t_theta[ i ] );
        _eq2[i] = MVCD__EQUIV::eq2( t_theta[ i ] );
        _eq3[i] = MVCD__EQUIV::eq3( t_theta[ i ] );
        _eq4[i] = MVCD__EQUIV::eq4( t_theta[ i ] );
        _eq5[i] = MVCD__EQUIV::eq5( t_theta[ i ] );
        _cos[i] = cos( t_theta[i] );
        dt_2[i] = ( pt[i] - eta ).sqrnorm();
        dt[i]   = sqrt( dt_2[i] );
        dt_3[i] = dt_2[i] * dt[i];
        dt_4[i] = dt_2[i] * dt_2[i];
        dt_5[i] = dt_4[i] * dt[i];
    }
    for( unsigned int i = 0 ; i < 3 ; ++i )
    {
        for( unsigned int j = 0 ; j < 3 ; ++j )
        {
            uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
            NN[i][j] = point_t::dot( N[i] , N[j] );
            JNN[i][j] = point_t::cross( pt[ (i+2)%3 ] - pt[ (i+1)%3 ] , N[j] );
        }
    }





    for( unsigned int t =  0 ; t < n_triangles ; ++t )
    {
        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        ///////           Check for a new test for Jacobians : d^2 < epsilon ??        ///////
        T determinant = point_t::dot( N[0] , pt[0]-eta );
        T sqrdist = determinant*determinant / (4 * point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] ).sqrnorm());
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
            point_is_on_triangle_support[t] = true;
        else
            point_is_on_triangle_support[t] = false;
        ///////////////////////////////////////////////////////////////////////////////////////

        if( point_is_on_triangle_support[t] )
        {
            // std::cout << "Special case needed for gradients computation" << std::endl;
            // GRADIENT SPECIAL CASE:
            {
                //////////////////                   CODE FOR GRADIENTS :                /////////////////////
                sum_gradW_on_triangles[t] = point_t(0,0,0);
                for( unsigned int i = 0 ; i < 3 ; ++i )
                {
                    point_t gradwTi( 0,0,0 );
                    for( unsigned int j = 0 ; j < 3 ; ++j )
                    {
                        gradwTi +=
                                (
                                    ( _eq2[j] * uu[i][j]               / (2 * dt[ (j+1)%3 ]   * dt[ (j+2)%3 ])   )
                                    + ( _eq1[j] * uu[j][j] * NN[i][j]  / (4 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ]) )
                                    + ( - NN[i][j]                     / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]) )
                                    ) * _nT_NT;
                    }
                    w_gradients[ cage_triangles[t][i] ] += gradwTi;
                    sum_gradW_on_triangles[t] += gradwTi;
                }
                //////////////////                   CODE FOR GRADIENTS                  /////////////////////
            }
        }
        else
        {
            // then compute Jm, compute B, and compute w_gradients:
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                _eq1[i] = MVCD__EQUIV::eq1( t_theta[ 3*t + i ] );
                _eq2[i] = MVCD__EQUIV::eq2( t_theta[ 3*t + i ] );
                _eq6[i] = MVCD__EQUIV::eq6( t_theta[ 3*t + i ] );
                _eq7[i] = MVCD__EQUIV::eq7( t_theta[ 3*t + i ] );
                _eq8[i] = MVCD__EQUIV::eq8( t_theta[ 3*t + i ] );
                _eq9[i] = MVCD__EQUIV::eq9( t_theta[ 3*t + i ] );
                _cos[i] = cos( t_theta[3*t + i] );
                dt_2[i] = ( pt[i] - eta ).sqrnorm();
                dt[i]   = sqrt( dt_2[i] );
                dt_3[i] = dt_2[i] * dt[i];
                dt_4[i] = dt_2[i] * dt_2[i];
                dt_5[i] = dt_4[i] * dt[i];
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                JN[i] = mat_t::vectorial( pt[(i+2)%3] - pt[(i+1)%3] );
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
                    NN[i][j] = point_t::dot( N[i] , N[j] );
                    JNN[i][j] = point_t::cross( pt[ (i+2)%3 ] - pt[ (i+1)%3 ] , N[j] );
                }
            }


            //////////////////                   CODE FOR GRADIENTS :                /////////////////////
            sum_gradW_on_triangles[t] = point_t(0,0,0);
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                point_t JmtNi(0,0,0);
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    JmtNi += -( _eq1[j] * NN[i][j] / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ]) ) * JNN[j][j]
                            - ( NN[i][j] / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]) ) * (2 * eta - pt[ (j+1)%3 ] - pt[ (j+2)%3 ])
                            - ( _eq2[j] / (2 * dt[ (j+1)%3 ] * dt[ (j+2)%3 ]) ) * JNN[j][i];
                }
                const point_t & gradwtiT = ( JmtNi + sum_W_on_triangles[t] * N[i] ) / determinant;
                w_gradients[ cage_triangles[t][i] ] += gradwtiT;
                sum_gradW_on_triangles[t] += gradwtiT;
            }
            //////////////////                   CODE FOR GRADIENTS                  /////////////////////
        }
    }




    for( unsigned int t =  0 ; t < n_triangles ; ++t )
    {
        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        ///////           Check for a new test for Hessians : d^3 < epsilon ??        ///////
        T determinant = point_t::dot( N[0] , pt[0]-eta );
        T sqrdist = determinant*determinant / (4 * point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] ).sqrnorm());
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
            point_is_on_triangle_support[t] = true;
        else
            point_is_on_triangle_support[t] = false;
        ///////////////////////////////////////////////////////////////////////////////////////

        if( point_is_on_triangle_support[t] )
        {
            // std::cout << "Special case needed for Hessians computation" << std::endl;
            //////////////////                   CODE FOR HESSIANS :                /////////////////////
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                point_t grad_dwTi(0,0,0);
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    grad_dwTi +=(
                            - _eq1[j] * uu[i][j] * JNN[j][j]                            / (2 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            - uu[i][j] * (2*eta-pt[(j+1)%3]-pt[(j+2)%3])                / (2 * dt_2[(j+1)%3]*dt_2[(j+2)%3])

                            + uu[j][j] * NN[i][j] * _eq4[j] * JNN[j][j]                 / (4 * dt_5[(j+1)%3]*dt_5[(j+2)%3])

                            + uu[j][j] * NN[i][j] * (2*eta-pt[(j+1)%3]-pt[(j+2)%3])     / (2 * dt_4[(j+1)%3]*dt_4[(j+2)%3])

                            - _eq1[j] * uu[j][j] * ( JNN[i][j] + JNN[j][i] )            / (4 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            - NN[i][j] * JNN[j][j]                                      / (1 * dt_4[(j+1)%3]*dt_4[(j+2)%3])

                            + _cos[j] * NN[i][j]*(2*eta-pt[(j+1)%3]-pt[(j+2)%3])        / (1 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            + ( JNN[i][j] + JNN[j][i] )                                 / (2 * dt_2[(j+1)%3]*dt_2[(j+2)%3])
                            );

                }
                w_Hessians[ cage_triangles[t][i] ] += ( mat_t::tensor( grad_dwTi , _nT_NT )  +  mat_t::tensor( _nT_NT , grad_dwTi ) );
            }
            //////////////////                   CODE FOR HESSIANS                  /////////////////////
        }
        else
        {
            // then compute Jm, compute B, and compute w_gradients:
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                _eq1[i] = MVCD__EQUIV::eq1( t_theta[ 3*t + i ] );
                _eq2[i] = MVCD__EQUIV::eq2( t_theta[ 3*t + i ] );
                _eq4[i] = MVCD__EQUIV::eq4( t_theta[ 3*t + i ] );
                _eq6[i] = MVCD__EQUIV::eq6( t_theta[ 3*t + i ] );
                _eq7[i] = MVCD__EQUIV::eq7( t_theta[ 3*t + i ] );
                _eq8[i] = MVCD__EQUIV::eq8( t_theta[ 3*t + i ] );
                _eq9[i] = MVCD__EQUIV::eq9( t_theta[ 3*t + i ] );
                _cos[i] = cos( t_theta[3*t + i] );
                dt_2[i] = ( pt[i] - eta ).sqrnorm();
                dt[i]   = sqrt( dt_2[i] );
                dt_3[i] = dt_2[i] * dt[i];
                dt_4[i] = dt_2[i] * dt_2[i];
                dt_5[i] = dt_4[i] * dt[i];
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                JN[i] = mat_t::vectorial( pt[(i+2)%3] - pt[(i+1)%3] );
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
                    NN[i][j] = point_t::dot( N[i] , N[j] );
                    JNN[i][j] = JN[i] * N[j];
                }
            }

            //////////////////                   CODE FOR HESSIANS :                /////////////////////
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                mat_t HwiT(0,0,0,0,0,0,0,0,0);
                determinant = point_t::dot( N[i] , pt[i]-eta );
                for( unsigned int c = 0 ; c < 3 ; ++c )
                {
                    point_t NitCTc(0,0,0);
                    point_t delta_c(0,0,0);
                    delta_c[c] = 1.0;
                    for( unsigned int j = 0 ; j < 3 ; ++j )
                    {
                        NitCTc +=
                                -_eq4[j] * JNN[j][j][c] * NN[i][j] * JNN[j][j]                      / (2 * dt_5[ (j+1)%3 ] * dt_5[ (j+2)%3 ])
                                -(2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * NN[i][j] * JNN[j][j]    / (1 * dt_4[ (j+1)%3 ] * dt_4[ (j+2)%3 ])
                                -_eq1[j] * JN[j] * ( mat_t::tensor( (JN[j].getCol(c)) , N[j] ) + mat_t::tensor( N[j] , (JN[j].getCol(c)) ) ) * N[i]
                                                                                                    / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                -JNN[j][j][c] * NN[i][j] * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] )    / (1 * dt_4[ (j+1)%3 ] * dt_4[ (j+2)%3 ])
                                +_cos[j] * NN[i][j] * (2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] )
                                                                                                    / (1 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                -( 2*NN[i][j] * delta_c + point_t::dot( (JN[j].getCol(c)) , N[i] ) * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] ) )
                                                                                                    / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ])
                                +_eq1[j] * JNN[j][j][c] * JNN[j][i]                                 / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                +(2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * JNN[j][i]               / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]);
                    }
                    HwiT.setRow( c , (NitCTc + sum_gradW_on_triangles[t][c]*N[i] + N[i][c]*sum_gradW_on_triangles[t] ) / determinant );
                }
                w_Hessians[ cage_triangles[t][i] ] += HwiT;
            }
            //////////////////                   CODE FOR HESSIANS                  /////////////////////
        }
    }

    point_t sumGradients(0,0,0);
    mat_t sumHessians(0,0,0,0,0,0,0,0,0);
    for(unsigned int v = 0; v < n_vertices; ++v)
    {
        sumGradients += w_gradients[v];
        sumHessians += w_Hessians[v];
    }

    gradients.resize(n_vertices);
    for(unsigned int v = 0; v < n_vertices; ++v)
        gradients[v] = w_gradients[v] / sumWeights - w_weights[v]*sumGradients/(sumWeights*sumWeights);

    Hessians.resize(n_vertices);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        Hessians[v] =
                w_Hessians[v] / sumWeights
                - w_weights[v] * sumHessians / (sumWeights*sumWeights)
                - ( mat_t::tensor( w_gradients[v] , sumGradients )
                    + mat_t::tensor( sumGradients , w_gradients[v] ) ) / (sumWeights*sumWeights)
                + 2 * w_weights[v] * mat_t::tensor(sumGradients,sumGradients) / (sumWeights*sumWeights*sumWeights);
    }
}











template< class point_t , class mat_t >
void computeCoordinatesAndDerivatives(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< point_t > & gradients ,
    std::vector< mat_t > & Hessians ,

    std::vector< point_t >  &  w_gradients,// unnormalized gradients
    std::vector< mat_t >  &  w_Hessians)// unnormalized Hessians
{
    typedef typename point_t::type_t    T;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();

    assert( cage_normals.size() == cage_triangles.size()   &&    "cage_normals.size() != cage_triangles.size()" );
    T epsilon = 0.000000001;

    weights.clear();
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[v] = 1.0;
            return;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    std::vector< T > w_weights(  n_vertices , 0.0 );

    unsigned int vid[3];
    T l[3];
    T theta[3] ;
    std::vector< T > t_theta( 3*n_triangles , 0.0 );
    T w[3];
    T c[3];
    T s[3];
    std::vector< bool > point_is_on_triangle_support(n_triangles, false);
    std::vector< T > sum_W_on_triangles(n_triangles , 0.0);


    for( unsigned int t = 0 ; t < n_triangles ; ++t )
    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            vid[i] =  cage_triangles[t][i];

        for( unsigned int i = 0 ; i <= 2 ; ++i )
            l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
            // maybe change the evaluation of theta, maybe keep sin_theta and express the Taylor expansion formula in terms of sin_theta
            t_theta[ 3*t + i ] = theta[i];
        }

        T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;

        if( M_PI - h < epsilon )
        {
            // eta is on the triangle t , use 2d barycentric coordinates :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
            {
                w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
            }
            sumWeights = w[0] + w[1] + w[2];

            weights.clear();
            weights.resize( n_vertices , 0.0 );
            weights[ vid[0] ] = w[0] / sumWeights;
            weights[ vid[1] ] = w[1] / sumWeights;
            weights[ vid[2] ] = w[2] / sumWeights;
            return;
        }

        // cos( psi_i ): corresponds to the dot product of the normals of the two adjacent spherical triangles.
        // they say in the paper that it is more robust to compute it that way...
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            c[ i ] = ( 2.0 * sin(h) * sin(h - theta[ i ]) ) / ( sin(theta[ (i+1) % 3 ]) * sin(theta[ (i+2) % 3 ]) ) - 1.0;

        T sign_Basis_u0u1u2 = 1;
        if( point_t::dot( point_t::cross(u[vid[0]] , u[vid[1]]) , u[vid[2]]) < 0.0 )
            sign_Basis_u0u1u2 = -1;

        // sin( psi_i ):
        for( unsigned int i = 0 ; i <= 2 ; ++i )
            s[ i ] = sign_Basis_u0u1u2 * sqrt( std::max( (T)(0.0) , (T)(1.0) - c[ i ] * c[ i ] ) );

        // test in original MVC paper: (they test if one angle psi is close to 0: it is "distance sensitive" in the sense that it does not
        // relate directly to the distance to the support plane of the triangle, and the farther away you go from the triangle, the worse it is)

        // simple test we suggest:
        // the determinant of the basis is 2*area(T)*d( eta , support(T) ), we can directly test for the distance to support plane of the triangle to be minimum
        T determinant = point_t::dot( cage_vertices[vid[0]] - eta , point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
        T sqrdist = determinant*determinant / (4 * point_t::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ).sqrnorm() );
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
        {
            point_is_on_triangle_support[t] = true;
            sum_W_on_triangles[t] = 0.0;
            continue;
        }
        // Note that epsilon should be multiplied by the diagonal of the cage model, in order to be correct
        // we cannot simply transform the model in the unit box: it does not change weights, but it scales Jacobians (linearly) and Hessians (quadratically)

        // maybe this should be reconsidered as well, as it does not increase at all the precision in the computations.
        // regarding the MVC original code: ( CAREFUL !!! THERE IS A DIRECT RENORMALIZATION HERE THAT HAS TO BE REMOVED !!!! )
        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            // original code:
            // w[ i ] = ( theta[ i ] - c[ (i+1)% 3 ]*theta[ (i+2) % 3 ] - c[ (i+2) % 3 ]*theta[ (i+1) % 3 ] ) / ( d[ vid[i] ] * sin( theta[ (i+1) % 3 ] ) * s[ (i+2) % 3 ] );
            // correct one (without renormalization)
            w[ i ] = ( theta[ i ] - c[ (i+1)% 3 ]*theta[ (i+2) % 3 ] - c[ (i+2) % 3 ]*theta[ (i+1) % 3 ] ) / ( 2 * d[ vid[i] ] * sin( theta[ (i+1) % 3 ] ) * s[ (i+2) % 3 ] );
        }

        sumWeights += ( w[0] + w[1] + w[2] );
        w_weights[ vid[0] ] += w[0];
        w_weights[ vid[1] ] += w[1];
        w_weights[ vid[2] ] += w[2];

        sum_W_on_triangles[t] = ( w[0] + w[1] + w[2] );
    }

    weights.resize(n_vertices,0.0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;

    // Computation of weights is done.



    std::vector< point_t > sum_gradW_on_triangles( n_triangles , point_t(0,0,0) );

    w_gradients.clear(  );
    w_gradients.resize( n_vertices , point_t(0,0,0) );
    w_Hessians.clear(  );
    w_Hessians.resize( n_vertices , mat_t(0,0,0,0,0,0,0,0,0) );

    // Code for Derivatives :
    point_t pt[3];  point_t N[3];
    T _eq1[3];      T _eq2[3];      T _eq3[3];      T _eq4[3];      T _eq5[3];      T _eq6[3];      T _eq7[3];      T _eq8[3];      T _eq9[3];      T _cos[3];
    T dt[3];        T dt_2[3];      T dt_3[3];      T dt_4[3];      T dt_5[3];
    T uu[3][3];
    T NN[3][3];
    point_t JNN[3][3];
    mat_t JN[3];




    const point_t  & NT = point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] );
    T NT_sqrnorm = point_t::dot( NT,NT );
    // ! NT_sqrnorm = 4 |T|^2 = NT_norm * 2 |T|
    const point_t & _nT_NT = - NT / (NT_sqrnorm);
    // _nT_NT = nT / (-2|T|)

    for( unsigned int i = 0 ; i < 3 ; ++i )
    {
        _eq1[i] = MVCD__EQUIV::eq1( t_theta[ i ] );
        _eq2[i] = MVCD__EQUIV::eq2( t_theta[ i ] );
        _eq3[i] = MVCD__EQUIV::eq3( t_theta[ i ] );
        _eq4[i] = MVCD__EQUIV::eq4( t_theta[ i ] );
        _eq5[i] = MVCD__EQUIV::eq5( t_theta[ i ] );
        _cos[i] = cos( t_theta[i] );
        dt_2[i] = ( pt[i] - eta ).sqrnorm();
        dt[i]   = sqrt( dt_2[i] );
        dt_3[i] = dt_2[i] * dt[i];
        dt_4[i] = dt_2[i] * dt_2[i];
        dt_5[i] = dt_4[i] * dt[i];
    }
    for( unsigned int i = 0 ; i < 3 ; ++i )
    {
        for( unsigned int j = 0 ; j < 3 ; ++j )
        {
            uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
            NN[i][j] = point_t::dot( N[i] , N[j] );
            JNN[i][j] = point_t::cross( pt[ (i+2)%3 ] - pt[ (i+1)%3 ] , N[j] );
        }
    }




    for( unsigned int t =  0 ; t < n_triangles ; ++t )
    {
        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        ///////           Check for a new test for Jacobians : d^2 < epsilon ??        ///////
        T determinant = point_t::dot( N[0] , pt[0]-eta );
        T sqrdist = determinant*determinant / (4 * point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] ).sqrnorm());
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
            point_is_on_triangle_support[t] = true;
        else
            point_is_on_triangle_support[t] = false;
        ///////////////////////////////////////////////////////////////////////////////////////

        if( point_is_on_triangle_support[t] )
        {
            // std::cout << "Special case needed for gradients computation" << std::endl;
            // GRADIENT SPECIAL CASE:
            {
                //////////////////                   CODE FOR GRADIENTS :                /////////////////////
                sum_gradW_on_triangles[t] = point_t(0,0,0);
                for( unsigned int i = 0 ; i < 3 ; ++i )
                {
                    point_t gradwTi( 0,0,0 );
                    for( unsigned int j = 0 ; j < 3 ; ++j )
                    {
                        gradwTi +=
                                (
                                    ( _eq2[j] * uu[i][j]               / (2 * dt[ (j+1)%3 ]   * dt[ (j+2)%3 ])   )
                                    + ( _eq1[j] * uu[j][j] * NN[i][j]  / (4 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ]) )
                                    + ( - NN[i][j]                     / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]) )
                                    ) * _nT_NT;
                    }
                    w_gradients[ cage_triangles[t][i] ] += gradwTi;
                    sum_gradW_on_triangles[t] += gradwTi;
                }
                //////////////////                   CODE FOR GRADIENTS                  /////////////////////
            }
        }
        else
        {
            // then compute Jm, compute B, and compute w_gradients:
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                _eq1[i] = MVCD__EQUIV::eq1( t_theta[ 3*t + i ] );
                _eq2[i] = MVCD__EQUIV::eq2( t_theta[ 3*t + i ] );
                _eq6[i] = MVCD__EQUIV::eq6( t_theta[ 3*t + i ] );
                _eq7[i] = MVCD__EQUIV::eq7( t_theta[ 3*t + i ] );
                _eq8[i] = MVCD__EQUIV::eq8( t_theta[ 3*t + i ] );
                _eq9[i] = MVCD__EQUIV::eq9( t_theta[ 3*t + i ] );
                _cos[i] = cos( t_theta[3*t + i] );
                dt_2[i] = ( pt[i] - eta ).sqrnorm();
                dt[i]   = sqrt( dt_2[i] );
                dt_3[i] = dt_2[i] * dt[i];
                dt_4[i] = dt_2[i] * dt_2[i];
                dt_5[i] = dt_4[i] * dt[i];
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                JN[i] = mat_t::vectorial( pt[(i+2)%3] - pt[(i+1)%3] );
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
                    NN[i][j] = point_t::dot( N[i] , N[j] );
                    JNN[i][j] = point_t::cross( pt[ (i+2)%3 ] - pt[ (i+1)%3 ] , N[j] );
                }
            }


            //////////////////                   CODE FOR GRADIENTS :                /////////////////////
            sum_gradW_on_triangles[t] = point_t(0,0,0);
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                point_t JmtNi(0,0,0);
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    JmtNi += -( _eq1[j] * NN[i][j] / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ]) ) * JNN[j][j]
                            - ( NN[i][j] / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]) ) * (2 * eta - pt[ (j+1)%3 ] - pt[ (j+2)%3 ])
                            - ( _eq2[j] / (2 * dt[ (j+1)%3 ] * dt[ (j+2)%3 ]) ) * JNN[j][i];
                }
                const point_t & gradwtiT = ( JmtNi + sum_W_on_triangles[t] * N[i] ) / determinant;
                w_gradients[ cage_triangles[t][i] ] += gradwtiT;
                sum_gradW_on_triangles[t] += gradwtiT;
            }
            //////////////////                   CODE FOR GRADIENTS                  /////////////////////
        }
    }




    for( unsigned int t =  0 ; t < n_triangles ; ++t )
    {
        for( unsigned int i = 0 ; i < 3 ; ++i )
            pt[i] = cage_vertices[ cage_triangles[t][i] ];
        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

        ///////           Check for a new test for Hessians : d^3 < epsilon ??        ///////
        T determinant = point_t::dot( N[0] , pt[0]-eta );
        T sqrdist = determinant*determinant / (4 * point_t::cross( pt[1] - pt[0] , pt[2] - pt[0] ).sqrnorm());
        T dist = sqrt( sqrdist );

        if( dist < epsilon )
            point_is_on_triangle_support[t] = true;
        else
            point_is_on_triangle_support[t] = false;
        ///////////////////////////////////////////////////////////////////////////////////////

        if( point_is_on_triangle_support[t] )
        {
            // std::cout << "Special case needed for Hessians computation" << std::endl;
            //////////////////                   CODE FOR HESSIANS :                /////////////////////
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                point_t grad_dwTi(0,0,0);
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    grad_dwTi +=(
                            - _eq1[j] * uu[i][j] * JNN[j][j]                            / (2 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            - uu[i][j] * (2*eta-pt[(j+1)%3]-pt[(j+2)%3])                / (2 * dt_2[(j+1)%3]*dt_2[(j+2)%3])

                            + uu[j][j] * NN[i][j] * _eq4[j] * JNN[j][j]                 / (4 * dt_5[(j+1)%3]*dt_5[(j+2)%3])

                            + uu[j][j] * NN[i][j] * (2*eta-pt[(j+1)%3]-pt[(j+2)%3])     / (2 * dt_4[(j+1)%3]*dt_4[(j+2)%3])

                            - _eq1[j] * uu[j][j] * ( JNN[i][j] + JNN[j][i] )            / (4 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            - NN[i][j] * JNN[j][j]                                      / (1 * dt_4[(j+1)%3]*dt_4[(j+2)%3])

                            + _cos[j] * NN[i][j]*(2*eta-pt[(j+1)%3]-pt[(j+2)%3])        / (1 * dt_3[(j+1)%3]*dt_3[(j+2)%3])

                            + ( JNN[i][j] + JNN[j][i] )                                 / (2 * dt_2[(j+1)%3]*dt_2[(j+2)%3])
                            );

                }
                w_Hessians[ cage_triangles[t][i] ] += ( mat_t::tensor( grad_dwTi , _nT_NT )  +  mat_t::tensor( _nT_NT , grad_dwTi ) );
            }
            //////////////////                   CODE FOR HESSIANS                  /////////////////////
        }
        else
        {
            // then compute Jm, compute B, and compute w_gradients:
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                _eq1[i] = MVCD__EQUIV::eq1( t_theta[ 3*t + i ] );
                _eq2[i] = MVCD__EQUIV::eq2( t_theta[ 3*t + i ] );
                _eq4[i] = MVCD__EQUIV::eq4( t_theta[ 3*t + i ] );
                _eq6[i] = MVCD__EQUIV::eq6( t_theta[ 3*t + i ] );
                _eq7[i] = MVCD__EQUIV::eq7( t_theta[ 3*t + i ] );
                _eq8[i] = MVCD__EQUIV::eq8( t_theta[ 3*t + i ] );
                _eq9[i] = MVCD__EQUIV::eq9( t_theta[ 3*t + i ] );
                _cos[i] = cos( t_theta[3*t + i] );
                dt_2[i] = ( pt[i] - eta ).sqrnorm();
                dt[i]   = sqrt( dt_2[i] );
                dt_3[i] = dt_2[i] * dt[i];
                dt_4[i] = dt_2[i] * dt_2[i];
                dt_5[i] = dt_4[i] * dt[i];
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                JN[i] = mat_t::vectorial( pt[(i+2)%3] - pt[(i+1)%3] );
            }
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                for( unsigned int j = 0 ; j < 3 ; ++j )
                {
                    uu[i][j] = point_t::dot( pt[ (i+2)%3 ] - pt[ (i+1)%3 ]  ,  pt[ (j+2)%3 ] - pt[ (j+1)%3 ] );
                    NN[i][j] = point_t::dot( N[i] , N[j] );
                    JNN[i][j] = JN[i] * N[j];
                }
            }

            //////////////////                   CODE FOR HESSIANS :                /////////////////////
            for( unsigned int i = 0 ; i < 3 ; ++i )
            {
                mat_t HwiT(0,0,0,0,0,0,0,0,0);
                determinant = point_t::dot( N[i] , pt[i]-eta );
                for( unsigned int c = 0 ; c < 3 ; ++c )
                {
                    point_t NitCTc(0,0,0);
                    point_t delta_c(0,0,0);
                    delta_c[c] = 1.0;
                    for( unsigned int j = 0 ; j < 3 ; ++j )
                    {
                        NitCTc +=
                                -_eq4[j] * JNN[j][j][c] * NN[i][j] * JNN[j][j]                      / (2 * dt_5[ (j+1)%3 ] * dt_5[ (j+2)%3 ])
                                -(2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * NN[i][j] * JNN[j][j]    / (1 * dt_4[ (j+1)%3 ] * dt_4[ (j+2)%3 ])
                                -_eq1[j] * JN[j] * ( mat_t::tensor( (JN[j].getCol(c)) , N[j] ) + mat_t::tensor( N[j] , (JN[j].getCol(c)) ) ) * N[i]
                                                                                                    / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                -JNN[j][j][c] * NN[i][j] * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] )    / (1 * dt_4[ (j+1)%3 ] * dt_4[ (j+2)%3 ])
                                +_cos[j] * NN[i][j] * (2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] )
                                                                                                    / (1 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                -( 2*NN[i][j] * delta_c + point_t::dot( (JN[j].getCol(c)) , N[i] ) * ( 2*eta - pt[(j+1)%3] - pt[(j+2)%3] ) )
                                                                                                    / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ])
                                +_eq1[j] * JNN[j][j][c] * JNN[j][i]                                 / (2 * dt_3[ (j+1)%3 ] * dt_3[ (j+2)%3 ])
                                +(2*eta[c]-pt[(j+1)%3][c]-pt[(j+2)%3][c]) * JNN[j][i]               / (2 * dt_2[ (j+1)%3 ] * dt_2[ (j+2)%3 ]);
                    }
                    HwiT.setRow( c , (NitCTc + sum_gradW_on_triangles[t][c]*N[i] + N[i][c]*sum_gradW_on_triangles[t] ) / determinant );
                }
                w_Hessians[ cage_triangles[t][i] ] += HwiT;
            }
            //////////////////                   CODE FOR HESSIANS                  /////////////////////
        }
    }

    point_t sumGradients(0,0,0);
    mat_t sumHessians(0,0,0,0,0,0,0,0,0);
    for(unsigned int v = 0; v < n_vertices; ++v)
    {
        sumGradients += w_gradients[v];
        sumHessians += w_Hessians[v];
    }

    gradients.resize(n_vertices);
    for(unsigned int v = 0; v < n_vertices; ++v)
        gradients[v] = w_gradients[v] / sumWeights - w_weights[v]*sumGradients/(sumWeights*sumWeights);

    Hessians.resize(n_vertices);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        Hessians[v] =
                w_Hessians[v] / sumWeights
                - w_weights[v] * sumHessians / (sumWeights*sumWeights)
                - ( mat_t::tensor( w_gradients[v] , sumGradients )
                    + mat_t::tensor( sumGradients , w_gradients[v] ) ) / (sumWeights*sumWeights)
                + 2 * w_weights[v] * mat_t::tensor(sumGradients,sumGradients) / (sumWeights*sumWeights*sumWeights);
    }
}






//<namespace DerivativesApproximations>
namespace DerivativesApproximations
{
//<namespace FiniteDifferences>
namespace FiniteDifferences
{
template< class point_t , class mat_t >
void computeCoordinatesAndDerivatives(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    typename point_t::type_t stencil_size ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< point_t > & gradients ,
    std::vector< mat_t > & Hessians)
{
    // Compute MVC at point eta
    // and we need to compute MVC on a 3x3x3 stencil centered in eta
    // (actually, we don't need to compute the values at the corners of the cube, so we need 19 evals in total)
    ::MVCoordinates::MVC3D::computeCoordinates(eta , cage_triangles , cage_vertices , cage_normals , weights );

    // [ +x , -x , +y , -y , +z , -z ]
    // [ +x+y , +x-y , -x+y , -x-y ; +x+z , +x-z , -x+z , -x-z ; +y+z , +y-z , -y+z , -y-z ]

    // [ +x , -x , +y , -y , +z , -z ] :
    std::vector< typename point_t::type_t > xweights , Xweights, yweights , Yweights, zweights , Zweights;
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,0) , cage_triangles , cage_vertices , cage_normals , xweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,0) , cage_triangles , cage_vertices , cage_normals , Xweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,0) , cage_triangles , cage_vertices , cage_normals , yweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,0) , cage_triangles , cage_vertices , cage_normals , Yweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,0,-1) , cage_triangles , cage_vertices , cage_normals , zweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,0,1) , cage_triangles , cage_vertices , cage_normals , Zweights );

    // [ +x+y , +x-y , -x+y , -x-y ; +x+z , +x-z , -x+z , -x-z ; +y+z , +y-z , -y+z , -y-z ] :
    std::vector< typename point_t::type_t >
            xyweights , Xyweights, xYweights , XYweights,
            yzweights , Yzweights, yZweights , YZweights,
            xzweights , xZweights, Xzweights , XZweights;
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,-1,0) , cage_triangles , cage_vertices , cage_normals , xyweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,-1,0)  , cage_triangles , cage_vertices , cage_normals , Xyweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,1,0)  , cage_triangles , cage_vertices , cage_normals , xYweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,1,0)   , cage_triangles , cage_vertices , cage_normals , XYweights );

    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,-1) , cage_triangles , cage_vertices , cage_normals , yzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,-1)  , cage_triangles , cage_vertices , cage_normals , Yzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,1)  , cage_triangles , cage_vertices , cage_normals , yZweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,1)   , cage_triangles , cage_vertices , cage_normals , YZweights );

    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,-1) , cage_triangles , cage_vertices , cage_normals , xzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,1)  , cage_triangles , cage_vertices , cage_normals , xZweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,-1)  , cage_triangles , cage_vertices , cage_normals , Xzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,1)   , cage_triangles , cage_vertices , cage_normals , XZweights );

    // apply Finite Differences scheme:
    gradients.resize( weights.size() );
    Hessians.resize( weights.size() );
    for( unsigned int w = 0 ; w < weights.size() ; ++w )
    {
        gradients[ w ][ 0 ] = ( Xweights[ w ] - xweights[ w ] ) / (2*stencil_size);
        gradients[ w ][ 1 ] = ( Yweights[ w ] - yweights[ w ] ) / (2*stencil_size);
        gradients[ w ][ 2 ] = ( Zweights[ w ] - zweights[ w ] ) / (2*stencil_size);

        Hessians[ w ](0 , 0) = ( Xweights[ w ] - 2*weights[ w ] + xweights[ w ] ) / (stencil_size*stencil_size);
        Hessians[ w ](1 , 1) = ( Yweights[ w ] - 2*weights[ w ] + yweights[ w ] ) / (stencil_size*stencil_size);
        Hessians[ w ](2 , 2) = ( Zweights[ w ] - 2*weights[ w ] + zweights[ w ] ) / (stencil_size*stencil_size);

        Hessians[ w ](0 , 1) = ( XYweights[ w ] - Xyweights[ w ] - xYweights[ w ] + xyweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](1 , 0) = Hessians[ w ](0 , 1);
        Hessians[ w ](0 , 2) = ( XZweights[ w ] - Xzweights[ w ] - xZweights[ w ] + xzweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](2 , 0) = Hessians[ w ](0 , 2);
        Hessians[ w ](1 , 2) = ( YZweights[ w ] - Yzweights[ w ] - yZweights[ w ] + yzweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](2 , 1) = Hessians[ w ](1 , 2);
    }
}


template< class point_t , class mat_t >
void computeCoordinatesAndDerivatives(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    typename point_t::type_t stencil_size ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< point_t > & gradients ,
    std::vector< mat_t > & Hessians ,
    std::vector< point_t > & w_gradients ,
    std::vector< mat_t > & w_Hessians)
{
    std::vector< typename point_t::type_t > w_weights;
    // Compute MVC at point eta
    // and we need to compute MVC on a 3x3x3 stencil centered in eta
    // (actually, we don't need to compute the values at the corners of the cube, so we need 19 evals in total)
    ::MVCoordinates::MVC3D::computeCoordinates(eta , cage_triangles , cage_vertices , cage_normals , weights , w_weights );

    // [ +x , -x , +y , -y , +z , -z ]
    // [ +x+y , +x-y , -x+y , -x-y ; +x+z , +x-z , -x+z , -x-z ; +y+z , +y-z , -y+z , -y-z ]

    // [ +x , -x , +y , -y , +z , -z ] :
    std::vector< typename point_t::type_t > xweights , Xweights, yweights , Yweights, zweights , Zweights;
    std::vector< typename point_t::type_t > w_xweights , w_Xweights, w_yweights , w_Yweights, w_zweights , w_Zweights;
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,0) , cage_triangles , cage_vertices , cage_normals , xweights , w_xweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,0) , cage_triangles , cage_vertices , cage_normals , Xweights , w_Xweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,0) , cage_triangles , cage_vertices , cage_normals , yweights , w_yweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,0) , cage_triangles , cage_vertices , cage_normals , Yweights , w_Yweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,0,-1) , cage_triangles , cage_vertices , cage_normals , zweights , w_zweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,0,1) , cage_triangles , cage_vertices , cage_normals , Zweights , w_Zweights );

    // [ +x+y , +x-y , -x+y , -x-y ; +x+z , +x-z , -x+z , -x-z ; +y+z , +y-z , -y+z , -y-z ] :
    std::vector< typename point_t::type_t >
            xyweights , Xyweights, xYweights , XYweights,
            yzweights , Yzweights, yZweights , YZweights,
            xzweights , xZweights, Xzweights , XZweights;
    std::vector< typename point_t::type_t >
            w_xyweights , w_Xyweights, w_xYweights , w_XYweights,
            w_yzweights , w_Yzweights, w_yZweights , w_YZweights,
            w_xzweights , w_xZweights, w_Xzweights , w_XZweights;
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,-1,0) , cage_triangles , cage_vertices , cage_normals , xyweights , w_xyweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,-1,0)  , cage_triangles , cage_vertices , cage_normals , Xyweights , w_Xyweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,1,0)  , cage_triangles , cage_vertices , cage_normals , xYweights , w_xYweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,1,0)   , cage_triangles , cage_vertices , cage_normals , XYweights , w_XYweights );

    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,-1) , cage_triangles , cage_vertices , cage_normals , yzweights , w_yzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,-1)  , cage_triangles , cage_vertices , cage_normals , Yzweights , w_Yzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,-1,1)  , cage_triangles , cage_vertices , cage_normals , yZweights , w_yZweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(0,1,1)   , cage_triangles , cage_vertices , cage_normals , YZweights , w_YZweights );

    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,-1) , cage_triangles , cage_vertices , cage_normals , xzweights , w_xzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(-1,0,1)  , cage_triangles , cage_vertices , cage_normals , xZweights , w_xZweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,-1)  , cage_triangles , cage_vertices , cage_normals , Xzweights , w_Xzweights );
    ::MVCoordinates::MVC3D::computeCoordinates(eta + stencil_size * point_t(1,0,1)   , cage_triangles , cage_vertices , cage_normals , XZweights , w_XZweights );

    // apply Finite Differences scheme:
    gradients.resize( weights.size() );
    Hessians.resize( weights.size() );
    for( unsigned int w = 0 ; w < weights.size() ; ++w )
    {
        gradients[ w ][ 0 ] = ( Xweights[ w ] - xweights[ w ] ) / (2*stencil_size);
        gradients[ w ][ 1 ] = ( Yweights[ w ] - yweights[ w ] ) / (2*stencil_size);
        gradients[ w ][ 2 ] = ( Zweights[ w ] - zweights[ w ] ) / (2*stencil_size);

        Hessians[ w ](0 , 0) = ( Xweights[ w ] - 2*weights[ w ] + xweights[ w ] ) / (stencil_size*stencil_size);
        Hessians[ w ](1 , 1) = ( Yweights[ w ] - 2*weights[ w ] + yweights[ w ] ) / (stencil_size*stencil_size);
        Hessians[ w ](2 , 2) = ( Zweights[ w ] - 2*weights[ w ] + zweights[ w ] ) / (stencil_size*stencil_size);

        Hessians[ w ](0 , 1) = ( XYweights[ w ] - Xyweights[ w ] - xYweights[ w ] + xyweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](1 , 0) = Hessians[ w ](0 , 1);
        Hessians[ w ](0 , 2) = ( XZweights[ w ] - Xzweights[ w ] - xZweights[ w ] + xzweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](2 , 0) = Hessians[ w ](0 , 2);
        Hessians[ w ](1 , 2) = ( YZweights[ w ] - Yzweights[ w ] - yZweights[ w ] + yzweights[ w ] ) / (4*stencil_size*stencil_size);
        Hessians[ w ](2 , 1) = Hessians[ w ](1 , 2);
    }

    w_gradients.resize( w_weights.size() );
    w_Hessians.resize( w_weights.size() );
    for( unsigned int w = 0 ; w < w_weights.size() ; ++w )
    {
        w_gradients[ w ][ 0 ] = ( w_Xweights[ w ] - w_xweights[ w ] ) / (2*stencil_size);
        w_gradients[ w ][ 1 ] = ( w_Yweights[ w ] - w_yweights[ w ] ) / (2*stencil_size);
        w_gradients[ w ][ 2 ] = ( w_Zweights[ w ] - w_zweights[ w ] ) / (2*stencil_size);

        w_Hessians[ w ](0 , 0) = ( w_Xweights[ w ] - 2*w_weights[ w ] + w_xweights[ w ] ) / (stencil_size*stencil_size);
        w_Hessians[ w ](1 , 1) = ( w_Yweights[ w ] - 2*w_weights[ w ] + w_yweights[ w ] ) / (stencil_size*stencil_size);
        w_Hessians[ w ](2 , 2) = ( w_Zweights[ w ] - 2*w_weights[ w ] + w_zweights[ w ] ) / (stencil_size*stencil_size);

        w_Hessians[ w ](0 , 1) = ( w_XYweights[ w ] - w_Xyweights[ w ] - w_xYweights[ w ] + w_xyweights[ w ] ) / (4*stencil_size*stencil_size);
        w_Hessians[ w ](1 , 0) = w_Hessians[ w ](0 , 1);
        w_Hessians[ w ](0 , 2) = ( w_XZweights[ w ] - w_Xzweights[ w ] - w_xZweights[ w ] + w_xzweights[ w ] ) / (4*stencil_size*stencil_size);
        w_Hessians[ w ](2 , 0) = w_Hessians[ w ](0 , 2);
        w_Hessians[ w ](1 , 2) = ( w_YZweights[ w ] - w_Yzweights[ w ] - w_yZweights[ w ] + w_yzweights[ w ] ) / (4*stencil_size*stencil_size);
        w_Hessians[ w ](2 , 1) = w_Hessians[ w ](1 , 2);
    }
}

}
//</namespace FiniteDifferences>









//<namespace TriCubic>
namespace TriCubic
{
// TODO: compute this stuff with more precision ?:
static
double M___TRICUBIC_333_CENTERED_POLYNOMIAL_INVERSE_ARRAY[64][64]=
{
  { -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 },
  { 0.0000813802083333148296162562473909929395 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333148296162562473909929395 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333356462979679690761258826 },
  { 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 },
  { -0.0000813802083333287074040640618477482349 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333330442127540038654842647 },
  { 0.0000813802083333148296162562473909929395 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333148296162562473909929395 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333425851918718763045035303 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 },
  { -0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111188209932265635870862752 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111188209932265635870862752 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111188209932265635870862752 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111084126523707027445198037 },
  { -0.0000813802083333287074040640618477482349 , 0.0000813802083333356462979679690761258826 , 0.0000813802083333287074040640618477482349 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333356462979679690761258826 , -0.0000813802083333356462979679690761258826 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333321768510160154619370587 },
  { 0.0000271267361111118820993226563587086275 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0000271267361111101473758466795516142156 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111101473758466795516142156 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , 0.0000271267361111101473758466795516142156 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111101473758466795516142156 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111101473758466795516142156 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111114484184536621569350245 },
  { 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 },
  { -0.0000813802083333356462979679690761258826 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333356462979679690761258826 , 0.0000813802083333356462979679690761258826 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333356462979679690761258826 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , -0.0000813802083333356462979679690761258826 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333356462979679690761258826 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333322852712332640123804595 },
  { -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000002168404344971008868015 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406249999998373696741271743348989 },
  { 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333332610531885009663710662 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332475006613448975656411 },
  { -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 },
  { 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0000813802083333356462979679690761258826 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , 0.0000813802083333356462979679690761258826 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111101473758466795516142156 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111117736791054078082652268 },
  { 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333330442127540038654842647 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333334914461501541360632928 },
  { -0.0000271267361111110147375846679551614216 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111110147375846679551614216 , 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332610531885009663710662 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333332610531885009663710662 , -0.0000271267361111112315780191650560482231 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111112315780191650560482231 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111112315780191650560482231 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000271267361111112315780191650560482231 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333333694734057495168144669 , -0.0000271267361111110689476932922303831219 },
  { 0.0000813802083333703407674875052180141211 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333148296162562473909929395 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333148296162562473909929395 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333425851918718763045035303 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333703407674875052180141211 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333148296162562473909929395 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333148296162562473909929395 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333425851918718763045035303 },
  { -0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111049432054187491303309798 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 },
  { -0.0000813802083333425851918718763045035303 , 0.0000813802083333425851918718763045035303 , 0.0000813802083333287074040640618477482349 , -0.0000813802083333356462979679690761258826 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333287074040640618477482349 , 0.0000813802083333287074040640618477482349 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333304421275400386548426468 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333356462979679690761258826 , -0.0000813802083333356462979679690761258826 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333339115744919922690314706 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333321768510160154619370587 },
  { 0.0000271267361111118820993226563587086275 , -0.0000813802083333287074040640618477482349 , 0.0000813802083333287074040640618477482349 , -0.0000271267361111101473758466795516142156 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000271267361111101473758466795516142156 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333321768510160154619370587 , -0.0000271267361111101473758466795516142156 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , 0.0000271267361111110147375846679551614216 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000271267361111101473758466795516142156 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111110147375846679551614216 },
  { -0.0000271267361111049432054187491303309798 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111049432054187491303309798 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111049432054187491303309798 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111118820993226563587086275 },
  { 0.0000090422453703775662603447926812805235 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703775662603447926812805235 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000090422453703775662603447926812805235 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703775662603447926812805235 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.1779785156250000000000000000000000000000 , 0.1779785156250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703775662603447926812805235 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703775662603447926812805235 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000090422453703775662603447926812805235 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703740968133928390670916997 },
  { 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111118820993226563587086275 },
  { -0.0000090422453703706273664408854529028758 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111110147375846679551614216 , 0.0000090422453703701936855718912511292729 },
  { -0.0000813802083333425851918718763045035303 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333287074040640618477482349 , 0.0000813802083333356462979679690761258826 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333287074040640618477482349 , -0.0000813802083333425851918718763045035303 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333287074040640618477482349 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333287074040640618477482349 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333356462979679690761258826 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333356462979679690761258826 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 },
  { 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111109063173674194047180208 },
  { 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333339115744919922690314706 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333330442127540038654842647 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333336947340574951681446692 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333331661854984084847330905 },
  { -0.0000271267361111118820993226563587086275 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111114484184536621569350245 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333339115744919922690314706 , -0.0000271267361111113399982364136064916238 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111114484184536621569350245 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111114484184536621569350245 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111112315780191650560482231 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333333694734057495168144669 , 0.0000271267361111111231578019165056048223 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111110147375846679551614216 },
  { 0.0000271267361111049432054187491303309798 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333287074040640618477482349 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0000271267361111084126523707027445198037 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111084126523707027445198037 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111101473758466795516142156 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111118820993226563587086275 },
  { -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703703021057891398015726736 },
  { -0.0000271267361111110147375846679551614216 , 0.0000271267361111110147375846679551614216 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111118820993226563587086275 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111110147375846679551614216 , -0.0000271267361111113399982364136064916238 , 0.0000271267361111114484184536621569350245 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111110147375846679551614216 , -0.0000271267361111112315780191650560482231 , 0.0000271267361111112315780191650560482231 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111110147375846679551614216 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111112451305463211248536481 },
  { 0.0000090422453703701936855718912511292729 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111110147375846679551614216 , -0.0000090422453703703021057891398015726736 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111112315780191650560482231 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111112315780191650560482231 , -0.0000090422453703701936855718912511292729 , 0.0000271267361111112315780191650560482231 , -0.0000271267361111112315780191650560482231 , 0.0000090422453703704105260063883520160744 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703704105260063883520160744 , 0.0000271267361111112315780191650560482231 , -0.0000271267361111112315780191650560482231 , 0.0000090422453703704105260063883520160744 , 0.0000271267361111112315780191650560482231 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111112315780191650560482231 , -0.0000271267361111112315780191650560482231 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333333694734057495168144669 , 0.0000271267361111111231578019165056048223 , 0.0000090422453703704105260063883520160744 , -0.0000271267361111112315780191650560482231 , 0.0000271267361111111231578019165056048223 , -0.0000090422453703703563158977640767943740 },
  { 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999861222121921855432447046 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999861222121921855432447046 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000034694469519536141888238 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000104083408558608425664715 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999973979147860347893583821 },
  { -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333356462979679690761258826 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333339115744919922690314706 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333335863138402466177012684 },
  { -0.0002441406250000017347234759768070944119 , 0.0002441406250000017347234759768070944119 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000004336808689942017736030 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000002168404344971008868015 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406249999982652765240231929055881 , 0.0002441406249999982652765240231929055881 , 0.0002441406250000004336808689942017736030 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000004336808689942017736030 , 0.0002441406250000004336808689942017736030 , 0.0002441406249999998915797827514495565993 , -0.0002441406249999998373696741271743348989 },
  { 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332610531885009663710662 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332610531885009663710662 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000004336808689942017736030 , 0.0002441406250000004336808689942017736030 , -0.0000813802083333334778936229980672578677 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333335321037316223424795680 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0000813802083333333559208785934480090418 },
  { -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333287074040640618477482349 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333425851918718763045035303 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333373810214439458832202945 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333356462979679690761258826 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 },
  { 0.0000271267361111049432054187491303309798 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111084126523707027445198037 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111084126523707027445198037 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111101473758466795516142156 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111113399982364136064916238 },
  { 0.0000813802083333339115744919922690314706 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333332610531885009663710662 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333321768510160154619370587 , -0.0000813802083333326105318850096637106617 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333326105318850096637106617 , -0.0000813802083333332610531885009663710662 , -0.0000813802083333329357925367553150408639 , 0.0000813802083333330442127540038654842647 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333332746057156570351764913 },
  { -0.0000271267361111110147375846679551614216 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111107978971501708542746201 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111112315780191650560482231 , 0.0000271267361111114484184536621569350245 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333330442127540038654842647 , -0.0000271267361111109063173674194047180208 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111114484184536621569350245 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333332610531885009663710662 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111118820993226563587086275 , 0.0000813802083333336947340574951681446692 , -0.0000813802083333332610531885009663710662 , 0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111114484184536621569350245 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333333152632971252415927665 , -0.0000271267361111110418426389800927722717 },
  { -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000030357660829594124152209 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406249999996747393482543486697978 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000004336808689942017736030 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000004336808689942017736030 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0002441406249999986989573930173946791911 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406249999998915797827514495565993 , 0.0002441406249999995663191310057982263970 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000004336808689942017736030 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0002441406249999995663191310057982263970 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406249999996205292396300734480974 },
  { 0.0000813802083333339115744919922690314706 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333340199947092408194748714 , -0.0000813802083333339115744919922690314706 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333330442127540038654842647 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333339115744919922690314706 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333335863138402466177012684 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333330442127540038654842647 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333329357925367553150408639 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333330442127540038654842647 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333331526329712524159276654 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333330442127540038654842647 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333330442127540038654842647 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333328273723195067645974632 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330848703354720719005400 },
  { 0.0002441406249999999457898913757247782996 , -0.0002441406250000001084202172485504434007 , -0.0002441406249999997831595655028991131985 , 0.0002441406250000000542101086242752217004 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000542101086242752217004 , -0.0002441406249999998915797827514495565993 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000542101086242752217004 , 0.0002441406250000000000000000000000000000 , -0.0002441406249999998915797827514495565993 , -0.0002441406250000002710505431213761085019 , 0.0002441406250000002710505431213761085019 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406249999999457898913757247782996 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000001084202172485504434007 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000542101086242752217004 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , 0.0002441406250000000542101086242752217004 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000001084202172485504434007 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000002168404344971008868015 , 0.0002441406249999998915797827514495565993 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , -0.0002441406249999998915797827514495565993 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0002441406249999998915797827514495565993 , -0.0002441406249999998915797827514495565993 , 0.0002441406249999998915797827514495565993 },
  { -0.0000813802083333332610531885009663710662 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333559208785934480090418 , 0.0000813802083333334236835143737920361673 , -0.0002441406250000000542101086242752217004 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333559208785934480090418 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333017107699691727873414 , -0.0000813802083333335049986773102048687178 , 0.0002441406250000002168404344971008868015 , -0.0002441406250000000542101086242752217004 , 0.0000813802083333332881582428131039819164 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333017107699691727873414 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333288158242813103981916 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333423683514373792036167 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333017107699691727873414 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333288158242813103981916 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333333152632971252415927665 , 0.0002441406249999999457898913757247782996 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333332746057156570351764913 , -0.0002441406249999999457898913757247782996 , 0.0002441406249999999457898913757247782996 , -0.0000813802083333332881582428131039819164 , 0.0000813802083333333288158242813103981916 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333332475006613448975656411 , 0.0002441406249999998915797827514495565993 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 },
  { 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333336947340574951681446692 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333336947340574951681446692 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0000813802083333339115744919922690314706 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333327189521022582141540624 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000004336808689942017736030 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0000813802083333342368351437379203616729 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330577652811599342896898 },
  { -0.0000271267361111123157801916505604822305 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111115568386709107073784253 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333335863138402466177012684 , 0.0000271267361111114484184536621569350245 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111114484184536621569350245 , 0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111107978971501708542746201 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000271267361111110689476932922303831219 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111109876325303558175505714 },
  { -0.0000813802083333334778936229980672578677 , 0.0000813802083333334236835143737920361673 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333332475006613448975656411 , 0.0002441406249999998915797827514495565993 , -0.0002441406250000000542101086242752217004 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000002710505431213761085019 , -0.0002441406249999995663191310057982263970 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000003794707603699265519026 , 0.0000813802083333330306602268477966788396 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333334643410958419984524426 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333333152632971252415927665 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333332746057156570351764913 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000542101086242752217004 , -0.0000813802083333333423683514373792036167 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333559208785934480090418 , -0.0000813802083333333423683514373792036167 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333332746057156570351764913 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000542101086242752217004 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333559208785934480090418 , 0.0000813802083333333017107699691727873414 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333331932905527206223439407 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333333017107699691727873414 , -0.0000813802083333331932905527206223439407 , 0.0002441406249999999457898913757247782996 , -0.0002441406249999999457898913757247782996 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999998915797827514495565993 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333334236835143737920361673 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333288158242813103981916 },
  { 0.0000271267361111111231578019165056048223 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333152632971252415927665 , -0.0000271267361111111096052747604367993972 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333332610531885009663710662 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332746057156570351764913 , -0.0000271267361111110553951661361615776968 , 0.0000813802083333332881582428131039819164 , -0.0000813802083333333017107699691727873414 , 0.0000271267361111110689476932922303831219 , -0.0000271267361111111231578019165056048223 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111111096052747604367993972 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000271267361111111028290111824023966847 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333423683514373792036167 , -0.0000271267361111111231578019165056048223 , -0.0000271267361111111231578019165056048223 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111111096052747604367993972 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000271267361111111028290111824023966847 , -0.0000813802083333333017107699691727873414 , 0.0000813802083333333017107699691727873414 , -0.0000271267361111110960527476043679939721 , 0.0000271267361111111028290111824023966847 , -0.0000813802083333333017107699691727873414 , 0.0000813802083333333017107699691727873414 , -0.0000271267361111110960527476043679939721 , -0.0000813802083333333288158242813103981916 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333423683514373792036167 , 0.0000813802083333333288158242813103981916 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333423683514373792036167 , -0.0000271267361111111096052747604367993972 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111111096052747604367993972 },
  { -0.0000813802083333425851918718763045035303 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333356462979679690761258826 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333495240857757835328811780 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333304421275400386548426468 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333304421275400386548426468 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333339115744919922690314706 },
  { 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000271267361111101473758466795516142156 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111109063173674194047180208 },
  { 0.0000813802083333321768510160154619370587 , -0.0000813802083333313094892780270583898528 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333330442127540038654842647 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333332610531885009663710662 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0000813802083333313094892780270583898528 , 0.0000813802083333336947340574951681446692 , -0.0000813802083333330442127540038654842647 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333336947340574951681446692 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333333830259329055856198920 },
  { -0.0000271267361111110147375846679551614216 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0000271267361111113399982364136064916238 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111110147375846679551614216 , 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333332610531885009663710662 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333331526329712524159276654 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333336947340574951681446692 , 0.0000813802083333336947340574951681446692 , -0.0000271267361111112315780191650560482231 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333332610531885009663710662 , -0.0000271267361111111502628562286432156725 },
  { 0.0000271267361111049432054187491303309798 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111084126523707027445198037 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , -0.0000813802083333287074040640618477482349 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , 0.0000813802083333287074040640618477482349 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333321768510160154619370587 , -0.0021972656250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , -0.0000271267361111101473758466795516142156 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111101473758466795516142156 , 0.0007324218750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111101473758466795516142156 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111110147375846679551614216 },
  { -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703697600047028970493556699 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , 0.0007324218750000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0197753906250000000000000000000000000000 , -0.0197753906250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111118820993226563587086275 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0065917968750000000000000000000000000000 , 0.0065917968750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703703021057891398015726736 },
  { -0.0000271267361111118820993226563587086275 , 0.0000271267361111110147375846679551614216 , 0.0000271267361111118820993226563587086275 , -0.0000271267361111114484184536621569350245 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111114484184536621569350245 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111112315780191650560482231 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333339115744919922690314706 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333339115744919922690314706 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111107978971501708542746201 , -0.0000271267361111113399982364136064916238 , 0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111112315780191650560482231 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111112451305463211248536481 },
  { 0.0000090422453703706273664408854529028758 , -0.0000271267361111110147375846679551614216 , 0.0000271267361111110147375846679551614216 , -0.0000090422453703704105260063883520160744 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0000090422453703706273664408854529028758 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111110147375846679551614216 , 0.0000090422453703704105260063883520160744 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111112315780191650560482231 , 0.0000271267361111110147375846679551614216 , -0.0000813802083333334778936229980672578677 , 0.0000813802083333334778936229980672578677 , -0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333334778936229980672578677 , 0.0000271267361111112315780191650560482231 , -0.0000090422453703701936855718912511292729 , 0.0000271267361111110147375846679551614216 , -0.0000271267361111110147375846679551614216 , 0.0000090422453703703021057891398015726736 , 0.0002441406250000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0000090422453703701936855718912511292729 , -0.0000271267361111110689476932922303831219 , 0.0000271267361111111231578019165056048223 , -0.0000090422453703703834209520762144052242 },
  { 0.0000813802083333321768510160154619370587 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333339115744919922690314706 , -0.0000813802083333321768510160154619370587 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333334778936229980672578677 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000813802083333334778936229980672578677 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333331526329712524159276654 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000004336808689942017736030 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406249999998915797827514495565993 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000013010426069826053208089 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406249999996747393482543486697978 , -0.0000813802083333326105318850096637106617 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333335863138402466177012684 , 0.0000813802083333326105318850096637106617 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333336947340574951681446692 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333343452553609864708050736 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333331661854984084847330905 },
  { -0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111114484184536621569350245 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111113399982364136064916238 , -0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333330442127540038654842647 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333331526329712524159276654 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000271267361111107978971501708542746201 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111110418426389800927722717 },
  { -0.0000813802083333335321037316223424795680 , 0.0000813802083333334778936229980672578677 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333333017107699691727873414 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333333152632971252415927665 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333332746057156570351764913 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333332746057156570351764913 , -0.0000813802083333331932905527206223439407 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406249999999457898913757247782996 , 0.0002441406250000000000000000000000000000 , -0.0002441406249999998373696741271743348989 , 0.0002441406249999998915797827514495565993 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000542101086242752217004 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406249999999457898913757247782996 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000001084202172485504434007 , 0.0002441406249999999457898913757247782996 , -0.0002441406250000000542101086242752217004 , 0.0000813802083333331932905527206223439407 , -0.0000813802083333333017107699691727873414 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333334101309872177232307422 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333333559208785934480090418 , -0.0000813802083333333423683514373792036167 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333559208785934480090418 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333334372360415298608415924 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333423683514373792036167 },
  { 0.0000271267361111111773679105407808265227 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333334236835143737920361673 , -0.0000271267361111111367103290725744102474 , -0.0000271267361111111231578019165056048223 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333333694734057495168144669 , 0.0000271267361111111367103290725744102474 , -0.0000271267361111111231578019165056048223 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111110960527476043679939721 , 0.0000271267361111110486189025581271749843 , -0.0000813802083333333017107699691727873414 , 0.0000813802083333332746057156570351764913 , -0.0000271267361111110825002204482991885470 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333017107699691727873414 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , -0.0000813802083333333288158242813103981916 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333423683514373792036167 , 0.0000813802083333332610531885009663710662 , -0.0002441406249999998915797827514495565993 , 0.0002441406249999999457898913757247782996 , -0.0000813802083333333017107699691727873414 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333965784600616544253171 , -0.0002441406250000001084202172485504434007 , 0.0002441406250000000542101086242752217004 , -0.0000813802083333333423683514373792036167 , -0.0000271267361111110825002204482991885470 , 0.0000813802083333332475006613448975656411 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111111096052747604367993972 , 0.0000271267361111111096052747604367993972 , -0.0000813802083333333559208785934480090418 , 0.0000813802083333333559208785934480090418 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111111028290111824023966847 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333288158242813103981916 , -0.0000271267361111111096052747604367993972 , -0.0000271267361111111231578019165056048223 , 0.0000813802083333334101309872177232307422 , -0.0000813802083333333423683514373792036167 , 0.0000271267361111110960527476043679939721 },
  { -0.0000271267361111118820993226563587086275 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111110147375846679551614216 , 0.0000813802083333339115744919922690314706 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0000271267361111114484184536621569350245 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111110147375846679551614216 , 0.0000813802083333339115744919922690314706 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333330442127540038654842647 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333336947340574951681446692 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333332610531885009663710662 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333330442127540038654842647 , 0.0002441406250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333336947340574951681446692 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , 0.0000271267361111118820993226563587086275 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000271267361111109063173674194047180208 , -0.0000813802083333339115744919922690314706 , 0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000813802083333331526329712524159276654 , 0.0000813802083333336947340574951681446692 , -0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000813802083333332610531885009663710662 , -0.0000271267361111113399982364136064916238 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000271267361111111638153833847120210976 },
  { 0.0000090422453703706273664408854529028758 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703705189462236369024594751 , -0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000090422453703701936855718912511292729 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703703021057891398015726736 , -0.0000271267361111110147375846679551614216 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111110147375846679551614216 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000813802083333334778936229980672578677 , 0.0021972656250000000000000000000000000000 , -0.0021972656250000000000000000000000000000 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333334778936229980672578677 , -0.0021972656250000000000000000000000000000 , 0.0021972656250000000000000000000000000000 , -0.0000813802083333333694734057495168144669 , -0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , -0.0000090422453703704105260063883520160744 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000090422453703704105260063883520160744 , 0.0000271267361111112315780191650560482231 , -0.0007324218750000000000000000000000000000 , 0.0007324218750000000000000000000000000000 , -0.0000271267361111111231578019165056048223 , -0.0000271267361111112315780191650560482231 , 0.0007324218750000000000000000000000000000 , -0.0007324218750000000000000000000000000000 , 0.0000271267361111111231578019165056048223 , 0.0000090422453703704105260063883520160744 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000090422453703703563158977640767943740 },
  { 0.0000271267361111110147375846679551614216 , -0.0000271267361111110689476932922303831219 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111110553951661361615776968 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333333694734057495168144669 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333332746057156570351764913 , 0.0000813802083333332610531885009663710662 , -0.0000813802083333333694734057495168144669 , -0.0000813802083333332610531885009663710662 , 0.0000813802083333332475006613448975656411 , -0.0000271267361111111231578019165056048223 , 0.0000271267361111110960527476043679939721 , 0.0000271267361111110689476932922303831219 , -0.0000271267361111110418426389800927722717 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333333152632971252415927665 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333332746057156570351764913 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000542101086242752217004 , 0.0000813802083333333423683514373792036167 , -0.0000813802083333333288158242813103981916 , -0.0000813802083333333559208785934480090418 , 0.0000813802083333333423683514373792036167 , 0.0000813802083333334778936229980672578677 , -0.0000813802083333333152632971252415927665 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333332475006613448975656411 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000542101086242752217004 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000001626303258728256651011 , -0.0000813802083333333559208785934480090418 , 0.0000813802083333333423683514373792036167 , 0.0000813802083333333559208785934480090418 , -0.0000813802083333334101309872177232307422 , -0.0000271267361111111773679105407808265227 , 0.0000271267361111110689476932922303831219 , 0.0000271267361111111163815383384712021098 , -0.0000271267361111110553951661361615776968 , 0.0000813802083333333017107699691727873414 , -0.0000813802083333333559208785934480090418 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333559208785934480090418 , -0.0000813802083333333017107699691727873414 , 0.0000813802083333333559208785934480090418 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333334372360415298608415924 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111111231578019165056048223 , -0.0000271267361111111163815383384712021098 , 0.0000271267361111111773679105407808265227 },
  { -0.0000090422453703703563158977640767943740 , 0.0000271267361111111773679105407808265227 , -0.0000271267361111111231578019165056048223 , 0.0000090422453703703427633706080079889489 , 0.0000271267361111110689476932922303831219 , -0.0000813802083333333694734057495168144669 , 0.0000813802083333333694734057495168144669 , -0.0000271267361111111231578019165056048223 , -0.0000271267361111110689476932922303831219 , 0.0000813802083333333694734057495168144669 , -0.0000813802083333334236835143737920361673 , 0.0000271267361111111638153833847120210976 , 0.0000090422453703703427633706080079889489 , -0.0000271267361111111502628562286432156725 , 0.0000271267361111111638153833847120210976 , -0.0000090422453703703834209520762144052242 , 0.0000271267361111111231578019165056048223 , -0.0000813802083333333152632971252415927665 , 0.0000813802083333333152632971252415927665 , -0.0000271267361111111096052747604367993972 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , -0.0000271267361111111028290111824023966847 , 0.0000813802083333333288158242813103981916 , -0.0000813802083333333423683514373792036167 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111110689476932922303831219 , 0.0000813802083333333152632971252415927665 , -0.0000813802083333333152632971252415927665 , 0.0000271267361111111096052747604367993972 , 0.0000813802083333333152632971252415927665 , -0.0002441406250000000000000000000000000000 , 0.0002441406250000000000000000000000000000 , -0.0000813802083333333288158242813103981916 , -0.0000813802083333333152632971252415927665 , 0.0002441406250000000000000000000000000000 , -0.0002441406250000000000000000000000000000 , 0.0000813802083333333288158242813103981916 , 0.0000271267361111110825002204482991885470 , -0.0000813802083333333423683514373792036167 , 0.0000813802083333333559208785934480090418 , -0.0000271267361111111231578019165056048223 , 0.0000090422453703703156583162958703780987 , -0.0000271267361111111502628562286432156725 , 0.0000271267361111111096052747604367993972 , -0.0000090422453703703563158977640767943740 , -0.0000271267361111110892764840263335912596 , 0.0000813802083333333559208785934480090418 , -0.0000813802083333333288158242813103981916 , 0.0000271267361111110960527476043679939721 , 0.0000271267361111111028290111824023966847 , -0.0000813802083333333288158242813103981916 , 0.0000813802083333333288158242813103981916 , -0.0000271267361111111096052747604367993972 , -0.0000090422453703703563158977640767943740 , 0.0000271267361111111231578019165056048223 , -0.0000271267361111111231578019165056048223 , 0.0000090422453703703766446884981800025116 }
};
/**************************** Computed with MATLAB ****************************/
/* source code: computeTriCubicCoeffs.m
function [A,invA] = computeTriCubicCoeffs(filename)
fid=fopen(filename, 'w');
A = zeros(64,64);
for(i = 0 : 3) for(j = 0 : 3) for(k = 0 : 3) for(x = 0 : 3) for(y = 0 : 3) for(z = 0 : 3)
  A( 1 + z + 4*y + 16*x , 1 + k + 4*j + 16*i ) = power(-3+2*x,i)*power(-3+2*y,j)*power(-3+2*z,k);
end end end end end end
invA = inv(A);
fprintf(fid,'{\n');
for(i=1:64)
    fprintf(fid,'  { ');
    for(j=1:64) fprintf(fid,'%3.40f',invA(i,j)); if(j<64) fprintf(fid,' , '); end end
    if(i<64) fprintf(fid,' },\n'); else fprintf(fid,' }\n'); end
end
fprintf(fid,'}');
fclose(fid);
*/


template< class point_t , class mat_t >
void computeCoordinatesAndDerivatives(
    point_t const & eta ,
    std::vector< std::vector< int > > const & cage_triangles ,
    std::vector< point_t > const & cage_vertices ,
    std::vector< point_t > const & cage_normals ,
    double stencil_size ,
    std::vector< typename point_t::type_t > & weights ,
    std::vector< point_t > & gradients ,
    std::vector< mat_t > & Hessians)
{
    // Compute MVC at point eta
    // and we need to compute MVC on a [-3stencil_size;3stencil_size]^3 stencil centered in eta
    // (actually, we don't need to compute the values at the corners of the cube, so we need 19 evals in total)
    ::MVCoordinates::MVC3D::computeCoordinates(eta , cage_triangles , cage_vertices , cage_normals , weights );

    // we suppose that we have a precomputed polynomial stencil [-3;3]^3 and that our point is in (0,0,0)
    std::vector< std::vector< typename point_t::type_t > > stencilWeights( 64 );

    for( int i = 0; i < 4 ; ++i )
        for( int j =0 ; j < 4 ; ++j )
            for( int k = 0 ; k < 4 ; ++k )
                ::MVCoordinates::MVC3D::computeCoordinates(
                        eta + stencil_size * point_t( -3+2*i , -3+2*j , -3+2*k )
                        , cage_triangles , cage_vertices , cage_normals ,
                        stencilWeights[ k + 4*j + 16*i ] );

    // CAREFUL WITH THE ORDERING !!!!!!!!

    gradients.resize( weights.size() );
    Hessians.resize( weights.size() );

    for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
    {
        //We process them coordinate per coordinate:
        double m__dummy_tricubic_a_coeffs[64];
        for( int i = 0; i < 4 ; ++i )
        {
            for( int j =0 ; j < 4 ; ++j )
            {
                for( int k = 0 ; k < 4 ; ++k )
                {
                    m__dummy_tricubic_a_coeffs[ k + 4*j + 16*i ] = 0.0;
                    for( int x = 0 ; x < 4 ; ++x )
                    {
                        for( int y = 0 ; y < 4 ; ++y )
                        {
                            for( int z = 0 ; z < 4 ; ++z )
                            {
                                m__dummy_tricubic_a_coeffs[ k + 4*j + 16*i ] +=
                                        M___TRICUBIC_333_CENTERED_POLYNOMIAL_INVERSE_ARRAY[ k + 4*j + 16*i ][ z + 4*y + 16*x ] *
                                        stencilWeights[ z + 4*y + 16*x ][v];
                            }
                        }
                    }
                }
            }
        }
        // At this point, m__dummy_tricubic_a_coeffs are the coeffs a(i,j,k) such that f(x,y,z) = \sum_{i,j,k}{ a(i,j,k)x^iy^jz^k
        // We get the derivatives in (0,0,0) easily
        gradients[v] = (1.0/stencil_size) *
                point_t(
                    m__dummy_tricubic_a_coeffs[ 0 + 4*0 + 16*1 ],//a_100 = f_x(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 0 + 4*1 + 16*0 ],//a_010 = f_y(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 1 + 4*0 + 16*0 ]//a_001 = f_z(0,0,0)
                    );

        Hessians[v] = (1.0/(stencil_size*stencil_size)) *
                mat_t(
                    2*m__dummy_tricubic_a_coeffs[ 0 + 4*0 + 16*2 ],//2a_200 = f_xx(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 0 + 4*1 + 16*1 ],//a_110 = f_xy(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 1 + 4*0 + 16*1 ],//a_101 = f_xz(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 0 + 4*1 + 16*1 ],//a_110 = f_yx(0,0,0)
                    2*m__dummy_tricubic_a_coeffs[ 0 + 4*2 + 16*0 ],//2a_020 = f_yy(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 1 + 4*1 + 16*0 ],//a_011 = f_yz(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 1 + 4*0 + 16*1 ],//a_101 = f_zx(0,0,0)
                    m__dummy_tricubic_a_coeffs[ 1 + 4*1 + 16*0 ],//a_011 = f_zy(0,0,0)
                    2*m__dummy_tricubic_a_coeffs[ 2 + 4*0 + 16*0 ]//2a_002 = f_zz(0,0,0)
                    );
    }
}
}
//</namespace TriCubic>
}
//</namespace DerivativesApproximations>



}
}









#endif // MVCOORDINATES3D_H
