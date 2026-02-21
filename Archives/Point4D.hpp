// Point4D.hpp

#ifndef POINT4D_HPP
#define POINT4D_HPP

//#define DEBUG

#include <iostream>
#include "Point.hpp"

typedef int coord_t; // this is the type of the coordinates of points
                     // in the rest of the code we will use coord_t
                     // instead of any particular type for the ease
                     // of future changes and portability

#define ROUND(a)       ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))
#define ROUND_L(a) ((a)>0 ? (unsigned long)((a)+0.5) : -(unsigned long)(0.5-(a)))

//**************************************************************
//
//	Class name: Point4D
//
//      4 dimensional version of Point class
//
//**************************************************************
template <class T = coord_t>
class Point4D
{
    friend class MesaPlot;
    friend class DrawGraph;
    
    friend ostream &operator<< <T>( ostream &output,
                                    const Point4D<T> &p );

    //multiplication by a number (made friend for commutativity)
    friend const Point4D<T>
    operator*<T>( const Point4D<T> &p, T constant );
    
    friend const Point4D<T>
    operator*<T>( T constant, const Point4D<T> &p);

    // scalar product with a vector with float coordinates
//      friend const float
//      operator*<T>( const Point4D<T> &p, Point<float> &fp );
    
//      friend const Point4D<T>
//      operator*<T>( float constant, const Point4D<T> &p);
    
    friend const Point4D<T>
    operator/<T>( const Point4D<T> &p, T constant );

    friend const Point4D<T>
    operator/<T>( const Point4D<T> &p, float constant );
    
public:
    //default costructor sets dim to 3 and all coordinates to 0
    Point4D();
    Point4D(T x, T y, T z, T w);
    Point4D(T *_coord, int _dim);
    
    const Point4D<T> &operator=(const Point4D<T> & rhs); // assign Point4D
    Point4D( const Point4D<T> & init );                  // copy constructor
    ~Point4D(){ delete [] coord; }                     // destructor
    

    // overloaded arithmetic operators
    const Point4D<T> operator+(const Point4D<T> &p) const;
    Point4D<T> operator+=(const Point4D<T> &p);
    const Point4D<T> operator-(const Point4D<T> &p) const;
    Point4D<T> operator-() const;
    Point4D<T> operator-=(const Point4D<T> &p);
    const T operator*(const Point4D<T> &p);            // scalar product
    const T operator*(const Point4D<T> &p) const;      // scalar product
    
    const Point4D<T> operator*=(T constant);
    const Point4D<T> operator*=(float constant);
    const Point4D<T> operator/=(T constant);           // pt <- pt/const
    const Point4D<T> operator/=(float constant);       // pt <- pt/const
    const Point4D<T> operator/=(double constant);       // pt <- pt/const
    
    bool operator<(const Point4D<T> &p) const;
    bool operator==( const Point4D<T> &right ) const;

    const unsigned long norm() const ; // the norm of the vector associated
                                // with the point
    float fnorm() const ; // the norm
    float fnorm2() const ; // the norm 
    unsigned long norm2() const;      // the square of the norm

    // scalar porduct of v-u and w-u, without the need to compute the differences
    // v-u, w-u first.
    const T scal_pr(const Point4D<T> &u,
                    const Point4D<T> &v,
                    const Point4D<T> &w);
    
   // access operators
    T getX() const { return coord[0]; }
    T getY() const { return coord[1]; }
    T getZ() const { return coord[2]; }
    T getW() const { return coord[3]; }

    // modification operators
    void setX(T x){ coord[0] = x; }
    void setY(T y){ coord[1] = y; }
    void setZ(T z){ coord[2] = z; }
    void setW(T w){ coord[3] = w; }    
    void set_to_zero(){coord[0] = 0;
                       coord[1] = 0;
                       coord[2] = 0;
                       coord[3] = 0;
    }
                   
private:
    size_t dim;  // dimension of the ambient space of the point
    T *coord;    // its coordinates 
};

//==============================================================
//
//          TEMPLATE CLASS MEMBER FUNCTION DEFINITIONS
//
//==============================================================

//**************************************************************
//
//      Point4D constructors
//
//**************************************************************
template<class T>
Point4D<T>::Point4D() :  dim(4), coord(new T[dim])
{
#ifdef DEBUG
    cout << "Enterting constructor for point of dim " << dim << endl;
#endif
    
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    coord[3] = 0;
        
#ifdef DEBUG
    cout << "Leaving constructor for point of dim " << dim << endl;
#endif 
}

template<class T>
Point4D<T>::Point4D(T x, T y, T z, T w) : dim(4), coord(new T[dim])
{
#ifdef DEBUG
    cout << "Enterting 3nd constructor for point of dim " << dim << endl;
#endif

    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
    coord[3] = w;
        
#ifdef DEBUG
    cout << "Leaving 3nd constructor for point of dim " << dim << endl;
#endif 
}

template<class T>
Point4D<T>::Point4D(T *_coord, int _dim) : dim(_dim), coord(new T[_dim])
{
#ifdef DEBUG
    cout << "Enterting 2nd constructor for point of dim " << dim << endl;
#endif
    
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = _coord[ind];
    
#ifdef DEBUG
    cout << "Leaving 2nd constructor for point of dim " << dim << endl;
#endif 
}

//**************************************************************
//
//	Method name : Copy constructor
//
//**************************************************************
template<class T>
Point4D<T>::Point4D( const Point4D<T> &init ) :
        dim(init.dim), coord(new T[dim])
{
#ifdef DEBUG
    cout << "enterting copy constructor for point of dim "
         << dim << endl;
#endif
   for ( size_t i = 0; i < dim; i++ )
       coord[ i ] = init.coord[ i ];  // copy init into object

#ifdef DEBUG
    cout << "Leaving copy constructor for Point4D of dim "
         << dim << endl;
#endif
}

//**************************************************************
//
//	Method name : assignment operator (=)
//
//**************************************************************
template<class T>
const Point4D<T> &Point4D<T>::operator=( const Point4D<T> &rhs )
        // const return avoids: ( a1 = a2 ) = a
{
    if ( &rhs != this ) // check for self-assignment
    {
        // for points of diffrent dimension, deallocate the original
        // lhs point, the allocate new lhs point.
        if( dim != rhs.dim )
        {
            delete [] coord; // reclaim space
            dim = rhs.dim;   // change the dim of the lhs to rhs.dim
            coord = new T[dim]; // create space for lhs point
        }
        for ( size_t i = 0; i < dim; i++ )
           coord[ i ] = rhs.coord[ i ];  // copy elements of rhs point
                                         // into lhs point
    }
    return *this;   // enables x = y = z;
}

//**************************************************************
//
//	Method name : <<
//
//	Description : Output Point4D with overloaded stream
//      insertion operator.
//
//**************************************************************
template<class T>
ostream &operator<<( ostream &output, const Point4D<T> &p )
{
    output << '(';
    for(size_t ind = 0; ind < p.dim - 1; ind++)
        output << p.coord[ind] << ", ";

    output << p.coord[p.dim-1] << ')';
    
    return output;   // enables cascaded calls
}

//**************************************************************
//
//	Method name : addition operator (+)
//
//**************************************************************
template<class T>
const Point4D<T> Point4D<T>::operator+(const Point4D<T> &p) const
{
//#define DEBUG
#ifdef DEBUG
    cout << "Enterting (+) for points  " << (*this)
         << " and " << p << endl;
#endif  
    T *new_coord = new T[dim];
    for(size_t ind = 0; ind < dim; ind++)
        new_coord[ind] = coord[ind] + p.coord[ind];

    Point4D<T> result(new_coord, dim);
    delete [] new_coord;

#ifdef DEBUG
    cout << "Leaving (+) for points  " << (*this)
         << " and " << p << endl;
#endif     
    return result;
}

//**************************************************************
//
//	Method name : addition assignment (+=)
//
//**************************************************************
template<class T>
Point4D<T> Point4D<T>::operator+=(const Point4D<T> &p)
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] += p.coord[ind];
    
    return *this;
}

//**************************************************************
//
//	Method name : subtraction operator (-)
//
//**************************************************************
template<class T>
const Point4D<T> Point4D<T>::operator-(const Point4D<T> &p) const
{
    T *new_coord = new T[dim];
    for(size_t ind = 0; ind < dim; ind++)
        new_coord[ind] = coord[ind] - p.coord[ind];

    Point4D result(new_coord, dim);
    delete [] new_coord;

    return result;
}
// negative
template<class T>
Point4D<T> Point4D<T>::operator-() const
{
    return (*this) * (-1);
}

//**************************************************************
//
//	Method name : subtraction assignment (-=)
//
//**************************************************************
template<class T>
Point4D<T> Point4D<T>::operator-=(const Point4D<T> &p)
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] -= p.coord[ind];
    
    return *this;
}

//**************************************************************
//
//	Method name : multiplication by numbers (*)
//
//**************************************************************
template<class T>
const Point4D<T> operator*( T constant, const Point4D<T> &p)
{
    T *new_coord = new T[p.dim];
    
    for(size_t ind = 0; ind < p.dim; ind++)
        new_coord[ind] = constant * p.coord[ind];

    Point4D<T> result(new_coord, p.dim);
    delete [] new_coord;

    return result;
}

//  template<class T>
//  const Point4D<T> operator*( float constant, const Point4D<T> &p)
//  {
//      T *new_coord = new T[p.dim];
    
//      for(size_t ind = 0; ind < p.dim; ind++)
//          new_coord[ind] = (T)(ROUND(constant * p.coord[ind]));

//      Point4D<T> result(new_coord, p.dim);
//      delete [] new_coord;

//      return result;
//  }

template<class T>
const Point4D<T> operator*( const Point4D<T> &p,
                              T constant)
{
    return constant * p;
}

//**************************************************************
//
//	Method name : (*=) 
//
//      Description : replace pt by pt * constant
//
//**************************************************************
template<class T>
const Point4D<T> Point4D<T>::operator*=( T constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] *= constant;
    
    return *this;
}

template<class T>
const Point4D<T> Point4D<T>::operator*=( float constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = (T)ROUND(constant * coord[ind]);
    
    return *this;
}


//**************************************************************
//
//	Method name : (/=) 
//
//      Description : replace pt by pt/constant
//
//**************************************************************
template<class T>
const Point4D<T> Point4D<T>::operator/=( T constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] /= constant;
    
    return *this;
}

template<class T>
const Point4D<T> Point4D<T>::operator/=( float constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = (T)((float)coord[ind]/constant);
    
    return *this;
}

template<class T>
const Point4D<T> Point4D<T>::operator/=( double constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = (T)(coord[ind]/constant);
    
    return *this;
}

//  template<class T>
//  const Point4D<T> Point4D<T>::operator/=( const float constant )
//  {
//      for(size_t ind = 0; ind < dim; ind++)
//          coord[ind] = (T)(coord[ind]/constant);
    
//      return *this;
//  }
//**************************************************************
//
//	Method name : (/) 
//
//      Description :  pt/constant
//
//**************************************************************
template<class T>
const Point4D<T> operator/( const Point4D<T> &p, T constant )
{
    T *new_coord = new T[p.dim];
    
    for(size_t ind = 0; ind < p.dim; ind++)
        new_coord[ind] = p.coord[ind] / constant;
    
    Point4D<T> result(new_coord, p.dim);
    delete [] new_coord;

    return result;
}

//  template<class T>
//  const Point4D<T> operator/( const Point4D<T> &p, float constant )
//  {
//      T *new_coord = new T[p.dim];
    
//      for(size_t ind = 0; ind < p.dim; ind++)
//          new_coord[ind] = (T)(p.coord[ind] / constant);
    
//      Point4D<T> result(new_coord, p.dim);
//      delete [] new_coord;

//      return result;
//  }

//**************************************************************
//
//	Method name : (*) scalar product 
//
//      Description : scalar product of vectors associated with
//      points
//
//**************************************************************
template<class T>
const T Point4D<T>::operator*(const Point4D<T> &p)
{
    T result = 0;
    for(size_t ind = 0; ind < p.dim; ind++)
        result += coord[ind] * p.coord[ind];

    return result;
}

template<class T>
const T Point4D<T>::operator*(const Point4D<T> &p) const
{
    T result = 0;
    for(size_t ind = 0; ind < p.dim; ind++)
        result += coord[ind] * p.coord[ind];

    return result;
}

//  template<class T>
//  const float operator*<T>( const Point4D<T> &p, Point<float> &fp )
//  {
//      float result = 0;
//      for(size_t ind = 0; ind < p.dim; ind++)
//          result += p.coord[ind] * fp.coord[ind];

//      return result;
//  }
//**************************************************************
//
//       scal_pr()
//       scalar porduct of uv = v-u and uw = w-u,
//       without the need to compute the differences v-u, w-u.
//
//**************************************************************
template<class T>
const T scal_pr(const Point4D<T> &u,
                const Point4D<T> &v,
                const Point4D<T> &w)
{
    T uv_x = v.getX() - u.getX();
    T uv_y = v.getY() - u.getY();
    T uv_z = v.getZ() - u.getZ();  

    T uw_x = w.getX() - u.getX();
    T uw_y = w.getY() - u.getY();
    T uw_z = w.getZ() - u.getZ();

    return uv_x * uw_x + uv_y * uw_y + uv_z * uw_z;
}


//**************************************************************
//
//	Method name : < operator
//
//**************************************************************
template<class T>
bool Point4D<T>::operator<(const Point4D<T> &p) const
{
    bool result = false;

    if( coord[0] < p.getX() ||
        (coord[0] == p.getX() && coord[1] < p.getY() ) ||
        (dim == 3 &&
         coord[0] == p.getX() &&
         coord[1] == p.getY() &&
         coord[2] < p.getZ()))
        result = true;

    return result;
}

//**************************************************************
//
//	Method name : norm
//
//	Description : the norm of a vector
//
//**************************************************************
template<class T>
const unsigned long Point4D<T>::norm() const
{
    return ROUND_L( sqrt( norm2() ) );
}

//**************************************************************
//
//	Method name : norm
//
//	Description : the norm of a vector
//
//**************************************************************
template<class T>
float Point4D<T>::fnorm() const
{
    float result = 0;
    for(size_t ind = 0; ind < dim; ind++)
        result += (float)coord[ind] * (float)coord[ind];

    return sqrt( result );
}

//**************************************************************
//
//	Method name : norm2
//
//	Description : the square of the norm
//
//**************************************************************
template<class T>
unsigned long Point4D<T>::norm2() const
{
    unsigned long result = 0;
    for(size_t ind = 0; ind < dim; ind++){
        result += (unsigned long)coord[ind] * (unsigned long)coord[ind];
//          debug("coord["<<ind<<"]="<<coord[ind] <<", coord["<<ind<<"]^2="<<coord[ind]*coord[ind] << ", result="<< result);
    }
    return result;
}

//**************************************************************
//
//	Method name : fnorm2
//
//	Description : the square of the norm
//
//**************************************************************
template<class T>
float Point4D<T>::fnorm2() const
{
    float result = 0;
    for(size_t ind = 0; ind < dim; ind++)
        result += (float)coord[ind] * (float)coord[ind];

    return result;
}

//**************************************************************
//

//	Method name : overloaded equality operator (==)
//
//**************************************************************
template<class T>
bool Point4D<T>::operator==( const Point4D<T> &right ) const
{
    bool result = false;
    size_t ind = 0;

    while( (result = ( coord[ind] == right.coord[ind])) && ind < dim )
    {
        ind++;
    }

    return result;
}

#endif
