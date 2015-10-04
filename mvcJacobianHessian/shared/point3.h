#ifndef POINT3_H__
#define POINT3_H__

/* file:                point3.h
 * description:         classes for 3D points and matrices
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



/*
  CAREFUL !!!! :
  p1 * p2 is commonly described as a dot product, and p1 % p2 is sometimes used as a cross product.
  don't use p1 * p2 as a dot product, use point3<T>::dot( p1 , p2 ) instead.
  don't use p1 % p2 as a cross product, use point3<T>::cross( p1 , p2 ) instead.
  Otherwise it won't give you the result you expect !
*/


#include <cassert>


/*
ABOUT SVD :
If you want to use the methods that project 2x2 or 3x3 matrices
onto the closest 2D/3D rotation, or similarity, or rotation
restricted to a 2D plane in 3D, you need to use a Singular Values Decomposition.
We provide implementation of this, using the GNU Scientific Library.
Uncomment #define COMP__USE_GSL_FOR_MAT33
*/
#define COMP__USE_GSL_FOR_MAT33

#ifdef COMP__USE_GSL_FOR_MAT33
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#endif


/*
ABOUT ARBITRARY PRECISION TYPES :
If you want to use mpfrc++ (c++ wrapper library for mpfr, which is the GNU arbitrary precision library) library,
you should define COMP__USE_MPFRCPLUSPLUS
You can obtain the source at http://www.holoborodko.com/pavel/mpfr/ for example.
In this project, arbitrary precision library was used to show empirical convergence of Finite Difference schemes towards
Mean Value Coordinates derivatives as introduced in the corresponding paper.
Note, that if you use mpfrc++ and define COMP__USE_GSL_FOR_MAT33 without defining COMP__USE_MPFRCPLUSPLUS, this source won't compile.
*/
//#define COMP__USE_MPFRCPLUSPLUS

#ifdef COMP__USE_MPFRCPLUSPLUS
#include <mpreal.h>
#endif






template< typename T >
class point3
{
public:
    typedef T               type_t;

    point3< T >( T x_ , T y_ , T z_) {v[0] = x_; v[1] = y_; v[2] = z_;}
    point3< T >(){v[0] = 0; v[1] = 0; v[2] = 0;}

    inline  T x() const {return v[0];}
    inline  T y() const {return v[1];}
    inline  T z() const {return v[2];}

    inline  T operator [] (unsigned int c) const
    {
        return v[c];
    }
    inline  T & operator [] (unsigned int c)
    {
        return v[c];
    }


    void operator = (const point3< T > & other)
    {
        v[0] = other.x();
        v[1] = other.y();
        v[2] = other.z();
    }
    void operator += (const point3< T > & other)
    {
        v[0] += other.x();
        v[1] += other.y();
        v[2] += other.z();
    }
    void operator -= (const point3< T > & other)
    {
        v[0] -= other.x();
        v[1] -= other.y();
        v[2] -= other.z();
    }

    // T2 is supposed to be scalar: int, float, double, mpreal (mpfrc++)
    template< class T2 >
    void operator *= (T2 s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }

    // T2 is supposed to be scalar: int, float, double, mpreal (mpfrc++)
    template< class T2 >
    void operator /= (T2 s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }


    void setZero()
    {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
    }

    point3<T> getOrthogonal() const
    {
        if( v[0] == 0 )
        {
            return point3<T>( 0 , v[2] , -v[1] );
        }
        else if( v[1] == 0 )
        {
            return point3<T>( v[2] , 0 , -v[0] );
        }

        return point3<T>( v[1] , -v[0] , 0 );
    }

    T sqrnorm() const
    {
        return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
    }
    T norm() const
    {
        return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }

    void normalize()
    {
        T _n = norm();
        v[0] /= _n;
        v[1] /= _n;
        v[2] /= _n;
    }

    inline static
    T dot( const point3< T > & p1 , const point3< T > & p2 )
    {
        return p1.x() * p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
    }

    inline static
    point3< T > cross( const point3< T > & p1 , const point3< T > & p2 )
    {
        return point3< T >(
                    p1.y() * p2.z() - p1.z() * p2.y(),
                    p1.z() * p2.x() - p1.x() * p2.z(),
                    p1.x() * p2.y() - p1.y() * p2.x()
                    );
    }
    inline static point3< T > Min ( const point3< T >  & p , const point3< T >  & p2 )
    {
        return point3< T >( std::min( p2[0] , p[0] ),
                            std::min( p2[1] , p[1] ),
                            std::min( p2[2] , p[2] ) );
    }

    inline static point3< T > Max ( point3< T > const & p , point3< T > const & p2 )
    {
        return point3< T >( std::max( p2[0] , p[0] ),
                            std::max( p2[1] , p[1] ),
                            std::max( p2[2] , p[2] ) );
    }

    inline static point3< T > Rand ( point3< T > const & p , point3< T > const & p2 )
    {
        T rx = p[0] + (p2[0] - p[0]) * (float)(rand()) / (float)( RAND_MAX );
        T ry = p[1] + (p2[1] - p[1]) * (float)(rand()) / (float)( RAND_MAX );
        T rz = p[2] + (p2[2] - p[2]) * (float)(rand()) / (float)( RAND_MAX );
        return point3<T>( rx , ry , rz );
    }

    inline static void glVertex( point3< T > const & p )
    {
        glVertex3f( p[0],p[1],p[2] );
    }

    inline static void glNormal( point3< T > const & p )
    {
        glNormal3f( p[0],p[1],p[2] );
    }

    inline static int getDimension()
    {
        return 3;
    }

private:
    T v[3];
};





// p1 - p2 :
template< typename T > inline point3< T > operator + (const point3< T > & p1 , const point3< T > & p2 )
{ return point3< T >( p1.x() + p2.x() , p1.y() + p2.y() , p1.z() + p2.z() ); }
template< typename T > inline point3< T > operator - (const point3< T > & p1 , const point3< T > & p2 )
{ return point3< T >( p1.x() - p2.x() , p1.y() - p2.y() , p1.z() - p2.z() ); }

// -p :
template< typename T > inline point3< T > operator - (const point3< T > & p2 )
{ return point3< T >( - p2.x() , - p2.y() , - p2.z() ); }

// p * scalar :
// T2 is supposed to be scalar: int, float, double, mpreal (mpfrc++)
template< typename T , typename T2 > inline point3< T > operator * (const point3< T > & p , T2 s)
{ return point3< T >( s*p.x() , s*p.y()  , s*p.z() ); }

// scalar * p :
// T2 is supposed to be scalar: int, float, double, mpreal (mpfrc++)
template< typename T , typename T2 > inline point3< T > operator * ( T2 s , const point3< T > & p )
{ return point3< T >( s*p.x() , s*p.y()  , s*p.z() ); }

// p / scalar :
// T2 is supposed to be scalar: int, float, double, mpreal (mpfrc++)
template< typename T , typename T2 > inline point3< T > operator / (const point3< T > & p , T2 s)
{ return point3< T >( p.x()/s , p.y()/s  , p.z()/s ); }

// cout :
template< typename T > inline std::ostream & operator << (std::ostream & s , point3< T > const & p)
{
    s << p[0] << " \t" << p[1] << " \t" << p[2];
    return s;
}

template< typename T >
inline double myAbs( const point3<T> & f )
{
    return f.norm();
}

typedef point3< float >    point3f;
typedef point3< double >   point3d;
typedef point3< long double >   point3ld;








template< typename T >
class mat33
{
private:
    T vals[9];
    // will be noted as :
    // 0 1 2
    // 3 4 5
    // 6 7 8

public:
    typedef T                       type_t;

    ////////////         CONSTRUCTORS          //////////////
    mat33<T>()
    {
        vals[0] = 0; vals[1] = 0; vals[2] = 0;
        vals[3] = 0; vals[4] = 0; vals[5] = 0;
        vals[6] = 0; vals[7] = 0; vals[8] = 0;
    }
    mat33<T>( T v1 , T v2 , T v3 , T v4 , T v5 , T v6 , T v7 , T v8 , T v9)
    {
        vals[0] = v1; vals[1] = v2; vals[2] = v3;
        vals[3] = v4; vals[4] = v5; vals[5] = v6;
        vals[6] = v7; vals[7] = v8; vals[8] = v9;
    }
    // provided if you need some "cast" from T2 to T :
    template< typename T2 > mat33<T>( const std::vector< T2 > & cc )
    {
        for(int i = 0 ; i < 9 ; ++i )
            vals[i]=cc[i];
    }

    void operator = (const mat33<T> & m)
    {
        for( unsigned int c = 0 ; c < 9 ; ++c )
            vals[c] = m(c);
    }

    void operator += (const mat33<T> & m)
    {
        for( unsigned int c = 0 ; c < 9 ; ++c )
            vals[c] += m(c);
    }

    // *= scalar :
    template< typename T2 >
    void operator *= (T2 s)
    {
        for( unsigned int c = 0 ; c < 9 ; ++c )
            vals[c] *= s;
    }

    // /= scalar :
    template< typename T2 >
    void operator /= (T2 s)
    {
        for( unsigned int c = 0 ; c < 9 ; ++c )
            vals[c] /= s;
    }

    /////////////         SET         ///////////////
    inline
    void setIdentity()
    {
        vals[0] = 1; vals[1] = 0; vals[2] = 0;
        vals[3] = 0; vals[4] = 1; vals[5] = 0;
        vals[6] = 0; vals[7] = 0; vals[8] = 1;
    }
    inline
    void setZero()
    {
        vals[0] = 0; vals[1] = 0; vals[2] = 0;
        vals[3] = 0; vals[4] = 0; vals[5] = 0;
        vals[6] = 0; vals[7] = 0; vals[8] = 0;
    }
    inline
    void set( const mat33<T> & other )
    {
        for(unsigned int i = 0 ; i < 9; ++i)
            vals[ i ] = other(i);
    }



    ////////        ACCESS TO COORDINATES      /////////
    T operator () (unsigned int i , unsigned int j) const
    { return vals[3*i + j]; }
    T & operator () (unsigned int i , unsigned int j)
    { return vals[3*i + j]; }
    T operator () (unsigned int v) const
    { return vals[v]; }
    T & operator () (unsigned int v)
    { return vals[v]; }

    T & getCoord(unsigned int i , unsigned int j)
    { return vals[3*i+j]; }
    T getCoord(unsigned int i , unsigned int j) const
    { return vals[3*i+j]; }
    T & getCoord(unsigned int i)
    { return vals[i]; }
    T getCoord(unsigned int i) const
    { return vals[i]; }





    ////////        ACCESS TO ROWS / COLUMNS      /////////
    point3<T> getRow( unsigned int i ) const
    {
        return point3<T>( vals[3*i] , vals[3*i+1] , vals[3*i+2] );
    }
    void setRow( unsigned int i , point3<T> p)
    {
        vals[3*i]   = p[0];
        vals[3*i+1] = p[1];
        vals[3*i+2] = p[2];
    }
    point3<T> getCol(unsigned int j) const
    {
        return point3<T>( vals[j] , vals[3+j] , vals[6+j] );
    }
    void setCol(unsigned int j , point3<T> p)
    {
        vals[j]   = p[0];
        vals[3+j] = p[1];
        vals[6+j] = p[2];
    }





    ////////        BASICS       /////////
    inline T sqrnorm()
    {
        return vals[0]*vals[0] + vals[1]*vals[1] + vals[2]*vals[2]
                + vals[3]*vals[3] + vals[4]*vals[4] + vals[5]*vals[5]
                + vals[6]*vals[6] +  vals[7]*vals[7] + vals[8]*vals[8];
    }

    inline T norm()
    { return std::sqrt( sqrnorm() ); }

    inline T determinant() const
    {
        return vals[0] * ( vals[4] * vals[8] - vals[7] * vals[5] )
                - vals[1] * ( vals[3] * vals[8] - vals[6] * vals[5] )
                + vals[2] * ( vals[3] * vals[7] - vals[6] * vals[4] );
    }

    inline T trace() const
    { return vals[0] + vals[4] + vals[8]; }




    ////////        TRANSPOSE       /////////
    inline
    void transpose()
    {
        T xy = vals[1] , xz = vals[2] , yz = vals[5];
        vals[1] = vals[3];
        vals[3] = xy;
        vals[2] = vals[6];
        vals[6] = xz;
        vals[5] = vals[7];
        vals[7] = yz;
    }
    mat33<T> getTranspose()
    {
        return mat33<T>(vals[0],vals[3],vals[6],vals[1],vals[5],vals[7],vals[2],vals[5],vals[8]);
    }



    inline
    mat33<T> getInverse() const
    {
        return mat33<T>
                (
                    vals[4]*vals[8] - vals[5]*vals[7], vals[2]*vals[7] - vals[1]*vals[8], vals[1]*vals[5] - vals[2]*vals[4],
                    vals[5]*vals[6] - vals[3]*vals[8], vals[0]*vals[8] - vals[2]*vals[6], vals[2]*vals[3] - vals[0]*vals[5],
                    vals[3]*vals[7] - vals[4]*vals[6], vals[1]*vals[6] - vals[0]*vals[7], vals[0]*vals[4] - vals[1]*vals[3]
                    ) / determinant();
    }



    // Symmetry
    inline
    T normalizedSymmetry()
    {
        return sqrt(
                    (
                        (vals[1]-vals[3])*(vals[1]-vals[3])+
                        (vals[2]-vals[6])*(vals[2]-vals[6])+
                        (vals[5]-vals[7])*(vals[5]-vals[7])
                        )/(
                        (vals[1])*(vals[1])+
                        (vals[2])*(vals[2])+
                        (vals[3])*(vals[3])+
                        (vals[5])*(vals[5])+
                        (vals[6])*(vals[6])+
                        (vals[7])*(vals[7])
                        ) );
    }



    ////////        HARMONICITY (not needed usually)       /////////
    void coutHessianHarmonicity()
    {
        std::cout << "xx + yy + zz = " << (getCoord(0,0) + getCoord(1,1) + getCoord(2,2)) << std::endl <<
                     "xy - yx : " << getCoord(0,1) - getCoord(1,0) << std::endl <<
                     "yz - zy : " << getCoord(1,2) - getCoord(2,1) << std::endl <<
                     "xz - zx : " << getCoord(0,2) - getCoord(2,0) << std::endl;
    }
    void enforceHessianHarmonicity()
    {
        vals[8] = -vals[0] - vals[4];
        vals[3] = vals[1];
        vals[5] = vals[7];
        vals[2] = vals[6];
    }
    inline
    T harmonicityError()
    {
        return std::abs( vals[8] + vals[0] + vals[4] ) + std::abs( vals[3] - vals[1] ) + std::abs( vals[5] - vals[7] ) + std::abs( vals[2] - vals[6] );
    }





    ////////        VECTOR PRODUCTS       /////////
    void setVectorProduct( const point3<T> & v1 , const point3<T> & v2 )
    {
        for( int i = 0 ; i < 3 ; ++i )
            for( int j = 0 ; j < 3 ; ++j )
                vals[ 3*i+j ] = v1[i]*v2[j];
    }

    void addVectorProduct( const point3<T> & v1 , const point3<T> & v2 )
    {
        for( int i = 0 ; i < 3 ; ++i )
            for( int j = 0 ; j < 3 ; ++j )
                vals[ 3*i+j ] += v1[i]*v2[j];
    }





    ////////////         ROTATION <-> AXIS/ANGLE         /////////////
    void getAxisAndAngleFromRotationMatrix( point3<T> & axis , T & angle )
    {
        angle = std::acos( (trace() - 1) / 2 );
        axis[0] = vals[7] - vals[5];
        axis[1] = vals[2] - vals[6];
        axis[2] = vals[3] - vals[1];
        axis.normalize();
    }

    inline static
    mat33<T> getRotationMatrixFromAxisAndAngle( const point3<T> & axis , T angle )
    {
        mat33<T> w = vectorial(axis);
        return Identity() + std::sin(angle) * w + (1 - std::cos(angle)) * w * w;
    }








    ////////////////////            STATIC STANDARD MATRICES          ////////////////////////////
    inline static mat33<T> Identity()
    {  return mat33<T>(1,0,0,0,1,0,0,0,1);  }

    inline static mat33<T> Zero()
    {  return mat33<T>(0,0,0,0,0,0,0,0,0);  }

    inline static mat33<T> Rand()
    {
        return mat33<T>( (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ),
                         (float)(rand()) / (float)( RAND_MAX ));
    }

    inline static mat33<T> Rand( double min_value , double max_value )
    {
        return mat33<T>( min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ),
                         min_value + (max_value - min_value) * ( (float)(rand()) / (float)( RAND_MAX ) ));
    }

    inline static mat33<T> RandRotation()
    {
        point3<T> axis(  (float)(rand()) / (float)( RAND_MAX )  ,
                         (float)(rand()) / (float)( RAND_MAX )  ,
                         (float)(rand()) / (float)( RAND_MAX )  );
        axis.normalize();
        T angle = 2.0 * M_PI * (float)(rand()) / (float)( RAND_MAX );

        return mat33<T>::getRotationMatrixFromAxisAndAngle( axis , angle );
    }

    // method defined as a standard Linear System Jacobian initialization (either Identity() or RandRotation())
    inline static
    mat33<T> LSJacobianInit()
    {
        return mat33<T>::Identity();
    //    return mat33<T>::RandRotation();
    }



#ifdef COMP__USE_GSL_FOR_MAT33
#ifdef COMP__USE_MPFRCPLUSPLUS
    // This is needed only because there is no automatic cast from mpreal to double types,
    // so, to some extent, we need to design a minimal set of functions to have compatibility with mpfrc++
    inline void gsl_matrix_set( gsl_matrix * u , unsigned int i , unsigned int j , mpfr::mpreal val ) const { gsl_matrix_set(u,i,j,val.toDouble()); }
#endif
#endif

    void SVD( mat33<T> & U , T & sx , T & sy , T & sz , mat33<T> & Vt ) const
    {
#ifdef COMP__USE_GSL_FOR_MAT33
        gsl_matrix * u = gsl_matrix_alloc(3,3);
        for(unsigned int i = 0 ; i < 3; ++i)
            for(unsigned int j = 0 ; j < 3; ++j)
                gsl_matrix_set( u , i , j , this->getCoord(i,j) );

        gsl_matrix * v = gsl_matrix_alloc(3,3);
        gsl_vector * s = gsl_vector_alloc(3);
        gsl_vector * work = gsl_vector_alloc(3);

        gsl_linalg_SV_decomp (u, v, s, work);

        sx = s->data[0];
        sy = s->data[1];
        sz = s->data[2];
        for(unsigned int i = 0 ; i < 3; ++i)
        {
            for(unsigned int j = 0 ; j < 3; ++j)
            {
                U(i,j) = gsl_matrix_get( u , i , j );
                Vt(i,j) = gsl_matrix_get( v , j , i );
            }
        }

        gsl_matrix_free(u);
        gsl_matrix_free(v);
        gsl_vector_free(s);
        gsl_vector_free(work);
#else
        assert( 0 && "You need to use a Linear Algebra Library in order to use mat33<T>::SVD !\nCurrent support: GNU Scientific Library only." );
#endif

        // Geometric interpretation of a Jacobian transformation matrix:
        // a transformation T is given as R.B.S.Bt, R = rotation , B = local basis (it's a rotation matrix), S = scales in the basis B
        // it can be obtained from the svd decomposition of T = U Sigma Vt :
        // B = V
        // S = Sigma
        // R = U.Vt
    }







    ///////////////////      Projections onto Rotations :     ////////////////////
    mat33<T> getRotationalPart()
    {
        mat33<T> U,Vt;
        T sx,sy,sz;
        SVD(U,sx,sy,sz,Vt);
        const mat33<T> & res = U*Vt;
        if( res.determinant() < 0 )
        {
            U(2) = -U(2);
            U(5) = -U(5);
            U(8) = -U(8);
            return U*Vt;
        }
        // else
        return res;
    }

    void setRotation()
    {
        mat33<T> U,Vt;
        T sx,sy,sz;
        SVD(U,sx,sy,sz,Vt);
        const mat33<T> & res = U*Vt;
        if( res.determinant() < 0 )
        {
            U(2) = -U(2);
            U(5) = -U(5);
            U(8) = -U(8);
            set(U*Vt);
            return;
        }
        // else
        set(res);
    }

    void setRotation( double weight )
    {
        mat33<T> U,Vt;
        T sx,sy,sz;
        SVD(U,sx,sy,sz,Vt);
        mat33<T> S( weight + (1-weight)*sx ,0,0,0,weight + (1-weight)*sy ,0,0,0, weight + (1-weight)*sz );
        const mat33<T> & res = U*S*Vt;
        if( res.determinant() < 0 )
        {
            S(8) = -S(8);
            set(U*S*Vt);
            return;
        }
        // else
        set(res);
    }

    void setSimilarity()
    {
        mat33<T> U,Vt;
        T sx,sy,sz;
        SVD(U,sx,sy,sz,Vt);
        mat33<T> S((sx+sy+sz)/3,0,0,0,(sx+sy+sz)/3,0,0,0,(sx+sy+sz)/3);
        const mat33<T> & res = U*S*Vt;
        if( res.determinant() < 0 )
        {
            S(8) = -S(8);
            set(U* S * Vt);
            return;
        }
        // else
        set(res);
    }

    // TO CHECK WHAT'S BEST
    // what is the L2 projection on the set of matrices of the form R . B2 . diag[1 1 x] . B2^T ?
    // is U . V^T . B2 . diag[1 1 x] . B2^T enough ?
    void setRotationOnTangentPlane( const mat33<T> & B2 )
    {   // B2 is assumed to be an orthogonal matrix with B2(2) = normal of the tangent plane
        // SVD :
        mat33<T> U,Vt;
        T sx,sy,sz;
        SVD(U,sx,sy,sz,Vt);

        if( (U*Vt).determinant() < 0 ){ U(2) = -U(2); U(5) = -U(5); U(8) = -U(8); }

        const point3<T> & VtB2_col2 = ( Vt * B2 ).getCol(2);
        T lambda = sx * VtB2_col2[0] * VtB2_col2[0]  +  sy * VtB2_col2[1] * VtB2_col2[1]  +  sz * VtB2_col2[2] * VtB2_col2[2];

        // if( (U*Vt).determinant() < 0 ){ lambda = -lambda; }

        mat33<T> NewScales = mat33<T>::Identity();
        NewScales(2,2) = lambda;

        set( U * Vt * B2 * transposeProduct01( NewScales , B2 ) );
    }

    void setRotationOnTangentPlane( const point3<T> & p_normal )
    {   // p_normal is assumed to be normalized
        // construction of B2 : 3x3 matrix that represents an orthogonal basis where B2(2) = p_normal
        point3<T> b1 = p_normal.getOrthogonal();
        b1.normalize();
        const point3<T> & b2 = point3<T>::cross( p_normal , b1 );

        mat33<T> B2;
        B2.setCol( 0 , b1 );
        B2.setCol( 1 , b2 );
        B2.setCol( 2 , p_normal );

        setRotationOnTangentPlane(B2);
    }






    //////////           Stupid products          /////////////
    // get the result of m1^T * m2:
    mat33<T> transposeProduct10( const mat33<T> & m1 , const mat33<T> & m2 )
    {
        return mat33<T>(
                    m1(0)*m2(0) + m1(3)*m2(3) + m1(6)*m2(6) ,
                    m1(0)*m2(1) + m1(3)*m2(4) + m1(6)*m2(7) ,
                    m1(0)*m2(2) + m1(3)*m2(5) + m1(6)*m2(8) ,
                    m1(1)*m2(0) + m1(4)*m2(3) + m1(7)*m2(6) ,
                    m1(1)*m2(1) + m1(4)*m2(4) + m1(7)*m2(7) ,
                    m1(1)*m2(2) + m1(4)*m2(5) + m1(7)*m2(8) ,
                    m1(2)*m2(0) + m1(5)*m2(3) + m1(8)*m2(6) ,
                    m1(2)*m2(1) + m1(5)*m2(4) + m1(8)*m2(7) ,
                    m1(2)*m2(2) + m1(5)*m2(5) + m1(8)*m2(8)
                    );
    }
    // get the result of m1 * m2^T:
    mat33<T> transposeProduct01( const mat33<T> & m1 , const mat33<T> & m2 )
    {
        return mat33<T>(
                    m1(0)*m2(0) + m1(1)*m2(1) + m1(2)*m2(2) ,
                    m1(0)*m2(3) + m1(1)*m2(4) + m1(2)*m2(5) ,
                    m1(0)*m2(6) + m1(1)*m2(7) + m1(2)*m2(8) ,
                    m1(3)*m2(0) + m1(4)*m2(1) + m1(5)*m2(2) ,
                    m1(3)*m2(3) + m1(4)*m2(4) + m1(5)*m2(5) ,
                    m1(3)*m2(6) + m1(4)*m2(7) + m1(5)*m2(8) ,
                    m1(6)*m2(0) + m1(7)*m2(1) + m1(8)*m2(2) ,
                    m1(6)*m2(3) + m1(7)*m2(4) + m1(8)*m2(5) ,
                    m1(6)*m2(6) + m1(7)*m2(7) + m1(8)*m2(8)
                    );
    }






    inline static
    mat33<T> tensor( const point3<T> & p1 , const point3<T> & p2 )
    {
        return mat33<T>(
                    p1.x()*p2.x() , p1.x()*p2.y() , p1.x()*p2.z(),
                    p1.y()*p2.x() , p1.y()*p2.y() , p1.y()*p2.z(),
                    p1.z()*p2.x() , p1.z()*p2.y() , p1.z()*p2.z());
    }

    inline static
    mat33<T> vectorial( const point3<T> & p )
    {
        return mat33<T>(
                    0       , -p.z()    , p.y()     ,
                    p.z()   , 0         , - p.x()   ,
                    - p.y() , p.x()     , 0
                    );
    }

};




template< class T >
mat33<T> operator + (const mat33<T> & m1 , const mat33<T> & m2)
{
    return mat33<T>( m1(0)+m2(0) , m1(1)+m2(1) , m1(2)+m2(2) , m1(3)+m2(3) , m1(4)+m2(4) , m1(5)+m2(5) , m1(6)+m2(6) , m1(7)+m2(7) , m1(8)+m2(8) );
}
template< class T >
mat33<T> operator - (const mat33<T> & m1 , const mat33<T> & m2)
{
    return mat33<T>( m1(0)-m2(0) , m1(1)-m2(1) , m1(2)-m2(2) , m1(3)-m2(3) , m1(4)-m2(4) , m1(5)-m2(5) , m1(6)-m2(6) , m1(7)-m2(7) , m1(8)-m2(8) );
}
template< class T >
mat33<T> operator - (const mat33<T> & m)
{
    return mat33<T>( -m(0) , -m(1) , -m(2) , -m(3) , -m(4) , -m(5) , -m(6) , -m(7) , -m(8) );
}





// CAREFUL WITH THOSE :
// at compilation, it cannot be used if T2 = point3<T> or mat33<T> , as they have to be understood as matrix-vector and matrix-matrix products.
//template< class T , class T2 >
//mat33<T> operator * (T2 s , const mat33<T> & m)
//{
//    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
//}
//template< class T , class T2 >
//mat33<T> operator * (const mat33<T> & m , T2 s)
//{
//    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
//}

template< class T , class T2 >
mat33<T> operator / (const mat33<T> & m , T2 s)
{
    return mat33<T>( m(0)/s , m(1)/s , m(2)/s , m(3)/s , m(4)/s , m(5)/s , m(6)/s , m(7)/s , m(8)/s );
}



template< class T >
mat33<T> operator * (int s , const mat33<T> & m)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (float s , const mat33<T> & m)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (double s , const mat33<T> & m)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (long double s , const mat33<T> & m)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
#ifdef COMP__USE_MPFRCPLUSPLUS
template< class T >
mat33<T> operator * (mpfr::mpreal s , const mat33<T> & m)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
#endif



template< class T >
mat33<T> operator * (const mat33<T> & m , int s)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (const mat33<T> & m , float s)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (const mat33<T> & m , double s)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
template< class T >
mat33<T> operator * (const mat33<T> & m , long double s)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
#ifdef COMP__USE_MPFRCPLUSPLUS
template< class T >
mat33<T> operator * (const mat33<T> & m , mpfr::mpreal s)
{
    return mat33<T>( m(0)*s , m(1)*s , m(2)*s , m(3)*s , m(4)*s , m(5)*s , m(6)*s , m(7)*s , m(8)*s );
}
#endif




template< typename T >
point3<T> operator * (const mat33<T> & m , const point3<T> & p) // computes m.p
{
    return point3<T>(
                m(0)*p[0] + m(1)*p[1] + m(2)*p[2],
                m(3)*p[0] + m(4)*p[1] + m(5)*p[2],
                m(6)*p[0] + m(7)*p[1] + m(8)*p[2]);
}
template< typename T >
point3<T> operator * (const point3<T> & p , const mat33<T> & m) // computes p^t . m = (m^t . p)^t
{
    return point3<T>(
                m(0)*p[0] + m(3)*p[1] + m(6)*p[2],
                m(1)*p[0] + m(4)*p[1] + m(7)*p[2],
                m(2)*p[0] + m(5)*p[1] + m(8)*p[2]);
}
template< typename T >
mat33<T> operator * (const mat33<T> & m1 , const mat33<T> & m2)
{
    return mat33<T>(
                m1(0)*m2(0) + m1(1)*m2(3) + m1(2)*m2(6) ,
                m1(0)*m2(1) + m1(1)*m2(4) + m1(2)*m2(7) ,
                m1(0)*m2(2) + m1(1)*m2(5) + m1(2)*m2(8) ,
                m1(3)*m2(0) + m1(4)*m2(3) + m1(5)*m2(6) ,
                m1(3)*m2(1) + m1(4)*m2(4) + m1(5)*m2(7) ,
                m1(3)*m2(2) + m1(4)*m2(5) + m1(5)*m2(8) ,
                m1(6)*m2(0) + m1(7)*m2(3) + m1(8)*m2(6) ,
                m1(6)*m2(1) + m1(7)*m2(4) + m1(8)*m2(7) ,
                m1(6)*m2(2) + m1(7)*m2(5) + m1(8)*m2(8)
                );
}




template< typename T > inline std::ostream & operator << (std::ostream & s , mat33< T > const & m)
{
    s << m(0) << " \t" << m(1) << " \t" << m(2) << std::endl << m(3) << " \t" << m(4) << " \t" << m(5) << std::endl << m(6) << " \t" << m(7) << " \t" << m(8) << std::endl;
    return s;
}



typedef mat33< float >  mat33f;
typedef mat33< double > mat33d;
typedef mat33< long double > mat33ld;



#endif
