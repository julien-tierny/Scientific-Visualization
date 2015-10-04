#ifndef MVCOORDINATES2D_H
#define MVCOORDINATES2D_H

/* file:                MVCoordinates2D.h
 * description:         functions computing the derivatives of the Mean Value Coordinates in 2D.
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
#include <cassert>


namespace MVCoordinates{
namespace MVC2D{

// MVC :
template< class T , class point_t >
bool computeCoordinates(
    const point_t & eta ,
    const std::vector< point_t > & cage_vertices ,
    const std::vector< std::pair< int,int > > & cage_oriented_edges ,
    std::vector< T > & weights
    )
{
    // The cage has to be a set of closed polylines, all oriented consistantly, preferably ClockWise (not sure what happens if not...)

    // If the point eta is on the cage, then we do not fill w_weights, as it is impossible to estimate gradients and Hessians around eta,
    // and therefore we do not need them.

    unsigned int n_vertices = cage_vertices.size() , n_edges = cage_oriented_edges.size();

    std::vector< T > w_weights( n_vertices , 0 );
    weights.clear();
    weights.resize( n_vertices , 0 );

    double epsilon = 0.01;
    double sqrEpsilon = epsilon * epsilon;

    for(unsigned int e = 0 ; e < n_edges; ++e)
    {
        // test if eta is on the edge, or on its support line
        int e0 = cage_oriented_edges[e].first, e1 = cage_oriented_edges[e].second;
        const point_t & p0 = cage_vertices[ e0 ];
        const point_t & p1 = cage_vertices[ e1 ];

        T dotProduct = (p0 - eta)*(p1 - eta);
        T sqrCosTheta = dotProduct * dotProduct / ( (p0-eta).sqrnorm() * (p1-eta).sqrnorm() );

        if( sqrCosTheta > 1.0 - sqrEpsilon ) // recall : cos(x)² ~ 1 - x²
        {
            if( dotProduct < 0 )
            {
                // consider it is on the edge // OK
                weights[ e0 ] = std::sqrt( (eta - p1).sqrnorm() / (p0 - p1).sqrnorm() );
                weights[ e1 ] = 1 - weights[ e0 ];
                w_weights.clear();
                // std::cout << "point ON THE EDGE " << e << " ->  LINEAR INTERPOLATION" << std::endl;
                return true;
            }
            else
            {
                // consider it is on the support line
                // the weights are unchanged
                // std::cout << "point in the support line of the edge " << e << " but not on it  ->  give 0 as weights" << std::endl;
            }
        }
        else
        {
            // general case
            const point_t & N0 = rotatePi_2( eta - p0 );
            const point_t & N1 = rotatePi_2( p1 - eta );
            const point_t & mE = ( N0 / N0.norm() + N1 / N1.norm() );

            T w0 = mE * N1 / ((p0 - eta)*N1);
            T w1 = mE * N0 / ((p1 - eta)*N0);
            w_weights[ e0 ] += w0;
            w_weights[ e1 ] += w1;
        }
    }

    T totalWeights = 0;
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        totalWeights += w_weights[v];

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v] = w_weights[v]/totalWeights;

    return false;
}



// MVC with structure that is needed for further computation of Jacobians and Hessians:

template< class T , class point_t >
void computeCoordinates(
    const point_t & eta ,
    const std::vector< point_t > & cage_vertices ,
    const std::vector< std::pair< int,int > > & cage_oriented_edges ,
    std::vector< T > & w_weights ,
    std::vector< std::pair<T,T> > & w_weights_on_edges ,
    std::vector< T > & weights
    )
{
    // The cage has to be a set of closed polylines, all oriented consistantly, preferably CCW (not sure what happen if not...)

    // If the point eta is on the cage, then we do not fill w_weights, as it is impossible to estimate gradients and Hessians around eta,
    // and therefore we do not need them.

    unsigned int n_vertices = cage_vertices.size() , n_edges = cage_oriented_edges.size();

    w_weights.clear();
    w_weights.resize( n_vertices , 0 );
    w_weights_on_edges.clear();
    w_weights_on_edges.resize(n_edges);
    weights.clear();
    weights.resize( n_vertices , 0 );

    double epsilon = 0.000001;
    double sqrEpsilon = epsilon * epsilon;

    for(unsigned int e = 0 ; e < n_edges; ++e)
    {
        // test if eta is on the edge, or on its support line
        int e0 = cage_oriented_edges[e].first, e1 = cage_oriented_edges[e].second;
        const point_t & p0 = cage_vertices[ e0 ];
        const point_t & p1 = cage_vertices[ e1 ];

        T dotProduct = (p0 - eta)*(p1 - eta);
        T sqrCosTheta = dotProduct * dotProduct / ( (p0-eta).sqrnorm() * (p1-eta).sqrnorm() );

        if( sqrCosTheta > 1.0 - sqrEpsilon ) // recall : cos(x)² ~ 1 - x²
        {
            if( dotProduct < 0 )
            {
                // consider it is on the edge
                weights[ e0 ] = std::sqrt( (eta - p1).sqrnorm() / (p0 - p1).sqrnorm() );
                weights[ e1 ] = 1 - weights[ e0 ];
                w_weights.clear();
                w_weights_on_edges.clear();
                return;
            }
            else
            {
                // consider it is on the support line
                // the weights are unchanged
                w_weights_on_edges[e] = std::pair<T,T>( 0,0 );
            }
        }
        else
        {
            // general case
            const point_t & N0 = rotatePi_2( eta - p0 );
            const point_t & N1 = rotatePi_2( p1 - eta );
            const point_t & mE = ( N0 / N0.norm() + N1 / N1.norm() );

            T w0 = mE * N1 / ((p0 - eta)*N1);
            T w1 = mE * N0 / ((p1 - eta)*N0);
            w_weights[ e0 ] += w0;
            w_weights[ e1 ] += w1;

            w_weights_on_edges[e] = std::pair<T,T>(w0,w1);
        }
    }

    T totalWeights = 0;
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        totalWeights += w_weights[v];

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v] = w_weights[v]/totalWeights;
}







template< class T , class point_t >
void computeGradients(
    const point_t & eta,
    const std::vector< point_t > & cage_vertices,
    const std::vector< std::pair<int,int> > & cage_oriented_edges,
    const std::vector< T > & w_weights,
    const std::vector< std::pair<T,T> > & w_weights_on_edges ,
    std::vector< point_t > & w_gradients,
    std::vector< std::pair<point_t,point_t> > & w_gradients_on_edges,
    std::vector< point_t > & gradients)
{
    unsigned int n_vertices = cage_vertices.size() , n_edges = cage_oriented_edges.size();

    w_gradients.clear();
    w_gradients.resize( n_vertices , point_t() );
    gradients.clear();
    gradients.resize( n_vertices , point_t() );

    w_gradients_on_edges.clear();
    w_gradients_on_edges.resize(n_edges);

    double epsilon = 0.000001;
    double sqrEpsilon = epsilon * epsilon;


    for(unsigned int e = 0 ; e < n_edges; ++e)
    {
        // test if eta is on the edge, or on its support line
        int e0 = cage_oriented_edges[e].first, e1 = cage_oriented_edges[e].second;
        const point_t & p0 = cage_vertices[ e0 ];
        const point_t & p1 = cage_vertices[ e1 ];

        T dotProduct = (p0 - eta)*(p1 - eta);
        T sqrCosTheta = dotProduct * dotProduct / ( (p0-eta).sqrnorm() * (p1-eta).sqrnorm() );

        if( sqrCosTheta > 1.0 - sqrEpsilon ) // recall : cos(x)² ~ 1 - x²
        {
            if( dotProduct < 0 )
            {
                // consider it is on the edge
                assert( 0 && "ERROR You should not compute Jacobians on the cage!!!" );
            }
            else
            {
                // consider it is on the support line
                const point_t & NE = rotatePi_2(p1 - p0);
                T sqrlengthE = NE.sqrnorm();
                const point_t & N0 = rotatePi_2( eta - p0 );
                const point_t & N1 = rotatePi_2( p1 - eta );
                T N0norm = N0.norm() , N1norm = N1.norm();

                const point_t & gradw0 = ( N0*N1/(2*N0norm*N0norm*N0norm) + 1/(2*N1norm)
                                           + 1/N0norm - 1/N1norm ) * NE / sqrlengthE;

                const point_t & gradw1 = ( 1/(2*N1norm) + N0*N1/(2*N1norm*N1norm*N1norm)
                                           - 1/N0norm + 1/N1norm ) * NE / sqrlengthE;

                w_gradients_on_edges[e] = std::pair<point_t,point_t>(gradw0,gradw1);

                w_gradients[e0] += gradw0;
                w_gradients[e1] += gradw1;
            }
        }
        else
        {
            // general case
            const point_t & N0 = rotatePi_2( eta - p0 );
            const point_t & N1 = rotatePi_2( p1 - eta );
            T N0norm = N0.norm() , N1norm = N1.norm();

            const point_t & RN0 = rotatePi_2(N0);
            const point_t & RN1 = rotatePi_2(N1);
            const point_t & BEtransposeN0 = - RN0/N0norm + RN0/N1norm
                    + ((eta - p0)*RN0) * (eta-p0) / (N0norm*N0norm*N0norm)
                    - ((eta - p1)*RN0) * (eta-p1) / (N1norm*N1norm*N1norm)
                    + (w_weights_on_edges[e].first + w_weights_on_edges[e].second) * N0;
            const point_t & BEtransposeN1 = - RN1/N0norm + RN1/N1norm
                    + ((eta - p0)*RN1) * (eta-p0) / (N0norm*N0norm*N0norm)
                    - ((eta - p1)*RN1) * (eta-p1) / (N1norm*N1norm*N1norm)
                    + (w_weights_on_edges[e].first + w_weights_on_edges[e].second) * N1;

            const point_t & gradw0 = BEtransposeN1 / ((p0-eta)*N1);
            const point_t & gradw1 = BEtransposeN0 / ((p1-eta)*N0);
            w_gradients_on_edges[e] = std::pair<point_t,point_t>(gradw0,gradw1);

            w_gradients[e0] += gradw0;
            w_gradients[e1] += gradw1;
        }
    }

    point_t sumGradients(0,0);
    T sumWeights = 0;
    for(unsigned int v = 0; v < n_vertices; ++v)
    {
        sumGradients += w_gradients[v];
        sumWeights += w_weights[v];
    }

    for(unsigned int v = 0; v < n_vertices; ++v)
        gradients[v] = w_gradients[v] / sumWeights - w_weights[v]*sumGradients/(sumWeights*sumWeights);
}







template< class T , class point_t , class mat_t >
void computeHessians(
    const point_t & eta,
    const std::vector< point_t > & cage_vertices,
    const std::vector< std::pair<int,int> > & cage_oriented_edges,
    const std::vector< T > & w_weights,
    const std::vector< std::pair<point_t,point_t> > & w_gradients_on_edges,
    const std::vector< point_t > & w_gradients,
    std::vector< mat_t > & w_Hessians,
    std::vector< mat_t > & Hessians)
{
    unsigned int n_vertices = cage_vertices.size() , n_edges = cage_oriented_edges.size();

    w_Hessians.clear();
    w_Hessians.resize( n_vertices , mat_t(0,0,0,0) );
    Hessians.clear();
    Hessians.resize( n_vertices );

    double epsilon = 0.000001;
    double sqrEpsilon = epsilon * epsilon;

    for(unsigned int e = 0 ; e < n_edges; ++e)
    {
        // test if eta is on the edge, or on its support line
        int e0 = cage_oriented_edges[e].first, e1 = cage_oriented_edges[e].second;
        const point_t & p0 = cage_vertices[ e0 ];
        const point_t & p1 = cage_vertices[ e1 ];

        T dotProduct = (p0 - eta)*(p1 - eta);
        T sqrCosTheta = dotProduct * dotProduct / ( (p0-eta).sqrnorm() * (p1-eta).sqrnorm() );

        if( sqrCosTheta > 1.0 - sqrEpsilon ) // recall : cos(x)² ~ 1 - x²
        {
            if( dotProduct < 0 )
            {
                // consider it is on the edge
                assert( 0 && "ERROR You should not compute Hessians on the cage!!!" );
            }
            else
            {
                // consider it is on the support line
                const point_t & NE = rotatePi_2(p1 - p0);
                T sqrlengthE = NE.sqrnorm();
                const point_t & N0 = rotatePi_2( eta - p0 );
                const point_t & N1 = rotatePi_2( p1 - eta );
                T N0norm = N0.norm() , N1norm = N1.norm();

                const point_t & RN0 = rotatePi_2(N0);
                const point_t & RN1 = rotatePi_2(N1);

                const point_t & graddw0 = (-3*RN0 - RN1)/(2 * N0norm*N0norm*N0norm)
                        +(-RN1)/(N1norm*N1norm*N1norm)
                        +(-3*(N0*N1) * RN0)/(2 * N0norm*N0norm*N0norm*N0norm*N0norm)
                        +(3 * RN1)/(2 * N1norm*N1norm*N1norm);

                const mat_t & Hw0 = mat_t::tensor(graddw0/sqrlengthE,NE);

                const point_t & graddw1 = (2*RN0)/(N0norm*N0norm*N0norm)
                        +(3*RN1 - RN0)/(2 * N1norm*N1norm*N1norm)
                        +(-3 * RN0)/(2 * N0norm*N0norm*N0norm)
                        +(3*(N0*N1)* RN1)/(2 * N1norm*N1norm*N1norm*N1norm*N1norm);

                const mat_t & Hw1 = mat_t::tensor(graddw1/sqrlengthE,NE);

                w_Hessians[e0] += Hw0;
                w_Hessians[e1] += Hw1;
            }
        }
        else
        {
            // general case: compute Cx Cy
            const point_t & deltax = point_t(1,0);
            const point_t & deltay = point_t(0,1);

            const point_t & p0eta = eta - p0;
            const point_t & p1eta = eta - p1;

            const point_t & N0 = rotatePi_2( eta - p0 );
            const point_t & N1 = rotatePi_2( p1 - eta );

            T p0etaNorm = p0eta.norm();
            T p0etaNorm3 = p0etaNorm*p0etaNorm*p0etaNorm;
            T p0etaNorm5 = p0etaNorm3*p0etaNorm*p0etaNorm;
            T p1etaNorm = p1eta.norm();
            T p1etaNorm3 = p1etaNorm*p1etaNorm*p1etaNorm;
            T p1etaNorm5 = p1etaNorm3*p1etaNorm*p1etaNorm;

            mat_t RPi_2(0,1,-1,0);
            mat_t Id(1,0,0,1);

            const mat_t & DxJm =  p1eta.x()*RPi_2 / p1etaNorm3
                    - p0eta.x()*RPi_2 / p0etaNorm3
                    - ( mat_t::tensor( deltax , p0eta ) + mat_t::tensor( p0eta , deltax ) ) / p0etaNorm3
                    + ( mat_t::tensor( deltax , p1eta ) + mat_t::tensor( p1eta , deltax ) ) / p1etaNorm3
                    + 3*p0eta.x()*( mat_t::tensor( p0eta , p0eta ) ) / p0etaNorm5
                    - 3*p1eta.x()*( mat_t::tensor( p1eta , p1eta ) ) / p1etaNorm5;

            const mat_t & DyJm =  p1eta.y()*RPi_2 / p1etaNorm3
                    - p0eta.y()*RPi_2 / p0etaNorm3
                    - ( mat_t::tensor( deltay , p0eta ) + mat_t::tensor( p0eta , deltay ) ) / p0etaNorm3
                    + ( mat_t::tensor( deltay , p1eta ) + mat_t::tensor( p1eta , deltay ) ) / p1etaNorm3
                    + 3*p0eta.y()*( mat_t::tensor( p0eta , p0eta ) ) / p0etaNorm5
                    - 3*p1eta.y()*( mat_t::tensor( p1eta , p1eta ) ) / p1etaNorm5;

            mat_t Cx = mat_t::tensor( deltax , w_gradients_on_edges[e].first + w_gradients_on_edges[e].second )
                    + DxJm + (w_gradients_on_edges[e].first.x() + w_gradients_on_edges[e].second.x()) * Id;
            mat_t Cy = mat_t::tensor( deltay , w_gradients_on_edges[e].first + w_gradients_on_edges[e].second )
                    + DyJm + (w_gradients_on_edges[e].first.y() + w_gradients_on_edges[e].second.y()) * Id;

            mat_t Hw0;
            Hw0.setRow(0, -(N1*Cx)/(p0eta*N1));
            Hw0.setRow(1, -(N1*Cy)/(p0eta*N1));

            mat_t Hw1;
            Hw1.setRow(0, -(N0*Cx)/(p1eta*N0));
            Hw1.setRow(1, -(N0*Cy)/(p1eta*N0));

            w_Hessians[e0] += Hw0;
            w_Hessians[e1] += Hw1;
        }
    }

    T sumWeights = 0;
    point_t sumGradients = point_t(0,0);
    mat_t sumHessians = mat_t(0,0,0,0);
    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        sumWeights += w_weights[v];
        sumGradients += w_gradients[v];
        sumHessians += w_Hessians[v];
    }

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        Hessians[v] = w_Hessians[v] / sumWeights - w_weights[v] * sumHessians / (sumWeights*sumWeights)
                - ( mat_t::tensor( w_gradients[v] , sumGradients ) + mat_t::tensor( sumGradients , w_gradients[v] ) ) / (sumWeights*sumWeights)
                + 2 * w_weights[v] * mat_t::tensor(sumGradients,sumGradients) / (sumWeights*sumWeights*sumWeights);
    }
}


template< class T , class point_t , class mat_t >
void buildCoordinatesAndDerivatives(
    const point_t & eta , // INPUT
    const std::vector< point_t > & cage_vertices , // INPUT
    const std::vector< std::pair<int,int> > & cage_edges , // INPUT
    std::vector< T > & weights , // OUTPUT
    std::vector< point_t > & gradients , // OUTPUT
    std::vector< mat_t > & Hessians) // OUTPUT
{
    /////////   NEEDED :   //////////
    std::vector< T > w_weights;
    std::vector< std::pair< T , T > > w_weights_on_edges;
    std::vector< mat22<T> > w_Hessians;
    std::vector< point2<T> > w_gradients;
    std::vector< std::pair< point2<T>,point2<T> > > w_gradients_on_edges;
    /////////////////////////////////
    MVCoordinates::MVC2D::computeCoordinates(eta, cage_vertices, cage_edges, w_weights, w_weights_on_edges, weights);
    MVCoordinates::MVC2D::computeGradients(eta, cage_vertices, cage_edges, w_weights, w_weights_on_edges, w_gradients, w_gradients_on_edges, gradients);
    MVCoordinates::MVC2D::computeHessians(eta, cage_vertices, cage_edges, w_weights, w_gradients_on_edges, w_gradients, w_Hessians, Hessians);

}

}
}



#endif // MVCOORDINATES2D_H

