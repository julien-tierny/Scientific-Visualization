#ifndef MVC_EQUIV_H
#define MVC_EQUIV_H


/* file:                MVC_Equiv.h
 * description:         safe evaluation of functions not defined in 0.
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
  In this file, we present the implementation of safe evaluation of many scalar functions, involved
  in the computation of the derivatives of the Mean Value Coordinates in the 3D case.
  See article for details.
*/


#define MVCD__EQUIV_TYPICAL_TOL         0.001




namespace MVCoordinates{
namespace MVCD__EQUIV{

// eq1(x) = (cos(x)*sin(x) - x) / sin(x)^3
template< class T >
T eq1( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =   - 2.0/3.0 - x2 / 5.0 - 17.0*x4/420.0 - 29.0*x6/4200.0
                - 1181.0*x8/1108800.0
                - 1393481.0*x10/9081072000.0
                - 763967.0*x12/36324288000.0
                - 133541.0*x14/48117888000.0;
//                - 3821869001.0*x16/10751460894720000.0
//                - 115665628927.0*x18/260185353652224000.0
//                - 8388993163723.0*x20/1538810520171724800000.0;
        return res;
    }
    T sx = sin(x);
    return (cos(x) * sx - x)/(sx*sx*sx);
}


// eq2(x) = x / sin(x)
template< class T >
T eq2( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =   1.0 + x2 / 6.0 + 7.0*x4/360.0 + 31.0*x6/15120.0
                + 127.0*x8/604800
                + 73.0*x10/3421440.0
                + 1414477.0*x12/653837184000.0
                + 8191.0*x14/37362124800.0;
//                + 16931177.0*x16/762187345920000.0
//                + 5749691557.0*x18/2554547108585472000.0
//                + 91546277357.0*x20/401428831349145600000.0;
        return res;
    }
    return x/sin(x);
}




// eq3(x) = (cos(x)-1) / sin(x)^2
template< class T >
T eq3( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =   -0.5 - x2 / 8.0 - x4/48.0 - 17.0*x6/5760.0
                - 31.0*x8/80640.0
                - 691.0*x10/14515200.0
                - 5461.0*x12/958003200.0
                - 929569.0*x14/1394852659200.0;
//                - 3202291.0*x16/41845579776000.0
//                - 221930581.0*x18/25609494822912000.0
//                - 4722116521.0*x20/4865804016353280000.0;
        return res;
    }
    T sx = sin(x);
    return (cos(x) - 1.0) / (sx*sx);
}




// eq4(x) = ( 2*cos(x)*sin(x)^3 + 3*(sin(x)*cos(x) - x) ) / sin(x)^5
template< class T >
T eq4( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =   -1.6 - 4.0*x2 / 7.0 - x4/7.0 - 211.0*x6/6930.0
                - 29509.0*x8/5045040.0
                - 157301.0*x10/151351200.0
                - 16079783.0*x12/92626934400.0
                - 61760113.0*x14/2239887686400.0;
//                - 30359523011.0*x16/7227370934784000.0
//                - 16640264468327.0*x18/26929184103005184000.0
//                - 617766523408427.0*x20/7001587866781347840000.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( 2.0*cx*sx*sx*sx + 3.0*(sx*cx - x) ) / (sx*sx*sx*sx*sx);
}






// eq5(x) = ( cos(x)*sin(x)^2*(1-2*cos(x)) - 2*cos(x)^2 + 2*cos(x) ) / sin(x)^4
template< class T >
T eq5( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
      //  T x16 = x14*x2;
      //  T x18 = x16*x2;
      //  T x20 = x18*x2;
        T res =  1.25 - x2 / 4.0 - 11.0*x4/192.0 - x6/90.0
                - 313.0*x8/161280.0
                - 2281.0*x10/7257600.0
                - 368107.0*x12/7664025600.0
                - 222247.0*x14/31701196800.0;
//                - 82516303.0*x16/83691159552000.0
//                - 1721992561.0*x18/12804747411456000.0
//                - 347979181721.0*x20/19463216065413120000.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( cx*sx*sx*(1.0-2.0*cx) - 2.0*cx*cx + 2.0*cx ) / (sx*sx*sx*sx);
}





// eq6(x) = ( 3*cos(x)*( cos(x)*x-sin(x) ) + cos(x)*sin(x)^3 ) / sin(x)^5
template< class T >
T eq6( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
      //  T x16 = x14*x2;
      //  T x18 = x16*x2;
      //  T x20 = x18*x2;
        T res =  -0.4 - x2 / 35.0 + 3.0*x4/140.0 + 1349.0*x6/138600.0
                + 267767.0*x8/100900800.0
                + 1752539.0*x10/3027024000.0
                + 204708709.0*x12/1852538688000.0
                + 862222247.0*x14/44797753728000.0;
//                + 1359125231539.0*x16/433642256087040000.0
//                + 260976933802873.0*x18/538583682060103680000.0
//                + 529743964972219.0*x20/7370092491348787200000.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( 3.0*cx*( cx*x-sx ) + cx*sx*sx*sx ) / (sx*sx*sx*sx*sx);
}







// eq7(x) = ( 3*( cos(x)*x-sin(x) ) + sin(x)^3 ) / sin(x)^3
template< class T >
T eq7( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =  - 2.0*x2 / 5.0 - 2.0*x4/21.0 - 4.0*x6/225.0
                - 2.0*x8/693.0
                - 2764.0*x10/6449625.0
                - 4.0*x12/66825.0
                - 28936.0*x14/3618239625.0;
//                - 87734.0*x16/84922212375.0
//                - 698444.0*x18/5373085843125.0
//                - 310732.0*x20/19405276970625.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( 3.0*( cx*x-sx ) + sx*sx*sx ) / (sx*sx*sx);
}




// eq8(x) = ( sin(x) - x*cos(x) ) * cos(x) / sin(x)^3
template< class T >
T eq8( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =  1.0 / 3.0 - x2 / 30.0 - 53.0*x4/2520.0 - 367.0*x6/75600.0
                - 5689.0*x8/6652800.0
                - 7198361.0*x10/54486432000.0
                - 1121539.0*x12/59439744000.0
                - 8117471.0*x14/3175780608000.0;
//                - 236480428279.0*x16/709596419051520000.0
//                - 5929710926423.0*x18/140500090972200960000.0
//                - 48228394603127.0*x20/9232863121030348800000.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( sx - x*cx )*cx / (sx*sx*sx);
}







// eq9(x) = ( sin(x) - x*cos(x) ) / sin(x)
template< class T >
T eq9( T x , T tol = MVCD__EQUIV_TYPICAL_TOL )
{
    if( abs(x) < tol )
    {
        T x2 = x*x;
        T x4 = x2*x2;
        T x6 = x4*x2;
        T x8 = x6*x2;
        T x10 = x8*x2;
        T x12 = x10*x2;
        T x14 = x12*x2;
     //   T x16 = x14*x2;
     //   T x18 = x16*x2;
     //   T x20 = x18*x2;
        T res =   x2 / 3.0 + x4/45.0 + 2.0*x6/945.0
                + x8/4725.0
                + 2.0*x10/93555.0
                + 1382.0*x12/638512875.0
                + 4.0*x14/18243225.0;
//                + 3617.0*x16/162820783125.0
//                + 87734.0*x18 / 38979295480125.0
//                + 349222.0*x20/1531329465290625.0;
        return res;
    }
    T sx = sin(x);
    T cx = cos(x);
    return ( sx - x*cx ) / (sx);
}

}
}




#endif // MVC_EQUIV_H
