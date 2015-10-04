#ifndef POINT2_H__
#define POINT2_H__

/* file:                point2.h
 * description:         classes for 2D points and matrices
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





#include <cassert>


#define __USE_GSL_FOR_MAT22

#ifdef __USE_GSL_FOR_MAT22
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#endif



template< typename T >
class point2
{
public:
    typedef T               type_t;

    point2< T >(){v[0] = 0; v[1] = 0;}
    point2< T >( T x_ , T y_ ){v[0] = x_; v[1] = y_;}

    T x() const {return v[0];}
    T y() const {return v[1];}

    T operator [] (unsigned int c) const
    {
        return v[c];
    }
    T & operator [] (unsigned int c)
    {
        return v[c];
    }

    void operator += (const point2< T > & other)
    {
        v[0] += other.x();
        v[1] += other.y();
    }
    void operator -= (const point2< T > & other)
    {
        v[0] -= other.x();
        v[1] -= other.y();
    }
    void operator *= (int s)
    {
        v[0] *= s;
        v[1] *= s;
    }
    void operator *= (float s)
    {
        v[0] *= s;
        v[1] *= s;
    }
    void operator *= (double s)
    {
        v[0] *= s;
        v[1] *= s;
    }
    void operator /= (int s)
    {
        v[0] /= s;
        v[1] /= s;
    }
    void operator /= (float s)
    {
        v[0] /= s;
        v[1] /= s;
    }
    void operator /= (double s)
    {
        v[0] /= s;
        v[1] /= s;
    }



    T dot( const point2< T > & p1 , const point2< T > & p2 )
    {
        return (p1.x()) * (p2.x()) + (p1.y()) * (p2.y());
    }

    T sqrnorm() const
    {
        return v[0]*v[0] + v[1]*v[1];
    }
    T norm() const
    {
        return std::sqrt(v[0]*v[0] + v[1]*v[1]);
    }

private:
    T v[2];
};


template< typename T >
inline point2< T > rotatePi_2( const point2< T > & p )// rotation of pi/2
{
    return point2< T >( -p.y() , p.x() );
}

template< typename T >
inline std::ostream & operator << (std::ostream & s , const point2< T > & p)
{
    s << p.x() << " " << p.y();
    return s;
}




template< typename T >
inline
point2< T > operator + (const point2< T > & p1 , const point2< T > & p2 )
{
    return point2< T >( p1.x() + p2.x() , p1.y() + p2.y() );
}





template< typename T >
inline
point2< T > operator - (const point2< T > & p1 , const point2< T > & p2 )
{
    return point2< T >( p1.x() - p2.x() , p1.y() - p2.y() );
}





template< typename T >
inline
point2< T > operator - (const point2< T > & p)
{
    return point2< T >( - p.x() , - p.y() );
}






template< typename T >
inline
point2< T > operator * (const point2< T > & p , int s)
{
    return point2< T >( s*p.x() , s*p.y() );
}
template< typename T >
inline
point2< T > operator * (const point2< T > & p , float s)
{
    return point2< T >( s*p.x() , s*p.y() );
}
template< typename T >
inline
point2< T > operator * (const point2< T > & p , double s)
{
    return point2< T >( s*p.x() , s*p.y() );
}






template< typename T >
inline
point2< T > operator / (const point2< T > & p , int s)
{
    return point2< T >( p.x()/s , p.y()/s );
}
template< typename T >
inline
point2< T > operator / (const point2< T > & p , float s)
{
    return point2< T >( p.x()/s , p.y()/s );
}
template< typename T >
inline
point2< T > operator / (const point2< T > & p , double s)
{
    return point2< T >( p.x()/s , p.y()/s );
}







template< typename T >
inline
point2< T > operator * ( int s , const point2< T > & p )
{
    return point2< T >( s*p.x() , s*p.y() );
}
template< typename T >
inline
point2< T > operator * ( float s , const point2< T > & p )
{
    return point2< T >( s*p.x() , s*p.y() );
}
template< typename T >
inline
point2< T > operator * ( double s , const point2< T > & p )
{
    return point2< T >( s*p.x() , s*p.y() );
}





template< typename T >
inline
T operator * (const point2< T > & p1 , const point2< T > & p2)
{
    return p1.x() * p2.x() + p1.y() * p2.y();
}





typedef point2< float >    point2f;
typedef point2< double >   point2d;







template< typename T >
class mat22
{
public:
    typedef T   type_t;

    mat22<T>()
    {
        vals[0] = 0;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
    }
    mat22<T>( T v1 , T v2 , T v3 , T v4 )
    {
        vals[0] = v1;
        vals[1] = v2;
        vals[2] = v3;
        vals[3] = v4;
    }

    void set( const mat22<T> & m )
    {
        vals[0] = m(0);
        vals[1] = m(1);
        vals[2] = m(2);
        vals[3] = m(3);
    }

    mat22<T> Id()
    {
        return mat22<T>(1,0,0,1);
    }
    mat22<T> Zero()
    {
        return mat22<T>(0,0,0,0);
    }
    void setZero()
    {
        vals[0] = 0;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
    }
    void setId()
    {
        vals[0] = 1;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 1;
    }
    void setRotation(T angle)
    {
        T cangle = std::cos(angle);
        T sangle = std::sin(angle);
        vals[0] = cangle;
        vals[1] = -sangle;
        vals[2] = sangle;
        vals[3] = cangle;
    }

    mat22<T> RPi_2()
    {
        return mat22<T>(0,-1,1,0);
    }

    void transpose()
    {
        T xy = vals[1];
        vals[1] = vals[2];
        vals[2] = xy;
    }
    mat22<T> getTranspose()
    {
        return mat22<T>(vals[0],vals[2],vals[1],vals[3]);
    }

    T operator () (unsigned int i , unsigned int j) const
    {
        return vals[2*i + j];
    }

    T operator () (unsigned int v) const
    {
        return vals[v];
    }
    T & operator () (unsigned int i , unsigned int j)
    {
        return vals[2*i + j];
    }

    T & operator () (unsigned int v)
    {
        return vals[v];
    }


    template< class point_t >
    point_t getRow( unsigned int i )
    {
        return point_t( vals[2*i] , vals[2*i+1] );
    }
    template< class point_t >
    void setRow( unsigned int i , const point_t & p)
    {
        vals[2*i]   = p.x();
        vals[2*i+1] = p.y();
    }

    template< class point_t >
    point_t getCol(unsigned int j)
    {
        return point_t( vals[j] , vals[2+j] );
    }
    template< class point_t >
    void setCol(unsigned int j , const point_t & p)
    {
        vals[j]   = p.x();
        vals[2+j] = p.y();
    }


    T determinant() const
    {
        return vals[0] * vals[3] - vals[1] * vals[2];
    }

#ifdef __USE_GSL_FOR_MAT22
    void SVD( mat22<T> & U , T & sx , T & sy , mat22<T> & V )
    {
        gsl_matrix * u = gsl_matrix_alloc(2,2);
        for(unsigned int i = 0 ; i < 4; ++i)
            u->data[i] = vals[i];
        gsl_matrix * v = gsl_matrix_alloc(2,2);
        gsl_vector * s = gsl_vector_alloc(2);
        gsl_vector * work = gsl_vector_alloc(2);

        gsl_linalg_SV_decomp (u,
                              v,
                              s,
                              work);

        sx = s->data[0];
        sy = s->data[1];
        U(0) = u->data[0];
        U(1) = u->data[1];
        U(2) = u->data[2];
        U(3) = u->data[3];
        V(0) = v->data[0];
        V(1) = v->data[1];
        V(2) = v->data[2];
        V(3) = v->data[3];

        gsl_matrix_free(u);
        gsl_matrix_free(v);
        gsl_vector_free(s);
        gsl_vector_free(work);

        // a transformation T is given as R.B.S.Bt, R = rotation , B = local basis (rotation matrix), S = scales in the basis B
        // it can be obtained from the svd decomposition of T = U Sigma Vt :
        // B = V
        // S = Sigma
        // R = U.Vt
    }

    mat22<T> getRotationalPart()
    {
        mat22<T> U,V;
        T sx,sy;
        SVD(U,sx,sy,V);
        V.transpose();
        const mat22<T> & res = U*V;
        if( res.determinant() < 0 )
        {
            U(0,0) = -U(0,0);
            U(0,1) = -U(0,1);
            return U*V;
        }
        // else
        return res;
    }

    void setRotation()
    {
        mat22<T> U,V;
        T sx,sy;
        SVD(U,sx,sy,V);
        V.transpose();
        const mat22<T> & res = U*V;
        if( res.determinant() < 0 )
        {
            U(0,0) = -U(0,0);
            U(0,1) = -U(0,1);
            set(U*V);
            return;
        }
        // else
        set(res);
    }

    void setSimilarity()
    {
        mat22<T> U,V;
        T sx,sy;
        SVD(U,sx,sy,V);
        mat22<T> S((sx+sy)/2,0,0,(sx+sy)/2);
        V.transpose();
        const mat22<T> & res = U*S*V;
        if( res.determinant() < 0 )
        {
            U(0,0) = -U(0,0);
            U(0,1) = -U(0,1);
            set(U* S * V);
            return;
        }
        // else
        set(res);
    }

#endif

    void operator += (const mat22<T> & m)
    {
        vals[0] += m(0);
        vals[1] += m(1);
        vals[2] += m(2);
        vals[3] += m(3);
    }
    void operator -= (const mat22<T> & m)
    {
        vals[0] -= m(0);
        vals[1] -= m(1);
        vals[2] -= m(2);
        vals[3] -= m(3);
    }

    T sqrnorm()
    {
        return vals[0]*vals[0] + vals[1]*vals[1] + vals[2]*vals[2] + vals[3]*vals[3];
    }

private:
    T vals[4];
    // will be noted as :
    // 0 1
    // 2 3
};

template< class T >
mat22<T> operator + (const mat22<T> & m1 , const mat22<T> & m2)
{
    return mat22<T>( m1(0)+m2(0) , m1(1)+m2(1) , m1(2)+m2(2) , m1(3)+m2(3) );
}
template< class T >
mat22<T> operator - (const mat22<T> & m1 , const mat22<T> & m2)
{
    return mat22<T>( m1(0)-m2(0) , m1(1)-m2(1) , m1(2)-m2(2) , m1(3)-m2(3) );
}

template< class T >
mat22<T> operator * (int s , const mat22<T> & m)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}
template< class T >
mat22<T> operator * (float s , const mat22<T> & m)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}
template< class T >
mat22<T> operator * (double s , const mat22<T> & m)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}


template< class T >
mat22<T> operator * (const mat22<T> & m , int s)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}
template< class T >
mat22<T> operator * (const mat22<T> & m , float s)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}
template< class T >
mat22<T> operator * (const mat22<T> & m , double s)
{
    return mat22<T>( s*m(0) , s*m(1) , s*m(2) , s*m(3) );
}


template< class T >
mat22<T> operator / (const mat22<T> & m , int s)
{
    return mat22<T>( m(0)/s , m(1)/s , m(2)/s , m(3)/s );
}
template< class T >
mat22<T> operator / (const mat22<T> & m , float s)
{
    return mat22<T>( m(0)/s , m(1)/s , m(2)/s , m(3)/s );
}
template< class T >
mat22<T> operator / (const mat22<T> & m , double s)
{
    return mat22<T>( m(0)/s , m(1)/s , m(2)/s , m(3)/s );
}



template< typename T >
point2<T> operator * (const mat22<T> & m , const point2<T> & p) // computes m.p
{
    return point2<T>(
                m(0,0)*p.x() + m(0,1)*p.y() ,
                m(1,0)*p.x() + m(1,1)*p.y() );
}
template< typename T >
point2<T> operator * (const point2<T> & p , const mat22<T> & m) // computes p^t . m = (m^t . p)^t
{
    return point2<T>(
                m(0,0)*p.x() + m(1,0)*p.y() ,
                m(0,1)*p.x() + m(1,1)*p.y() );
}

template< typename T >
mat22<T> operator * (const mat22<T> & m1 , const mat22<T> & m2)
{
    return mat22<T>( m1(0,0)*m2(0,0) + m1(0,1)*m2(1,0) , m1(0,0)*m2(0,1) + m1(0,1)*m2(1,1),
                m1(1,0)*m2(0,0) + m1(1,1)*m2(1,0) , m1(1,0)*m2(0,1) + m1(1,1)*m2(1,1));
}

template< typename T >
mat22<T> tensor( const point2<T> & p1 , const point2<T> & p2 )
{
    return mat22<T>(p1.x()*p2.x() , p1.x()*p2.y(),
               p1.y()*p2.x() , p1.y()*p2.y());
}

typedef mat22< float >  mat22f;
typedef mat22< double > mat22d;





#endif
