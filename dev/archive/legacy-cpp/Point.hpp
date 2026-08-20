// Point.hpp

 #ifndef POINT_HPP
 #define POINT_HPP

//#define DEBUG_POINT

#include <iostream>

typedef int coord_t; // this is the type of the coordinates of points
                     // in the rest of the code we will use coord_t
                     // instead of any particular type for the ease
                     // of future changes and portability

#define ROUND(a)       ((a)>0 ? (int)((a)+0.5) : -(int)(0.5-(a)))
#define ROUND_L(a) ((a)>0 ? (unsigned long)((a)+0.5) : -(unsigned long)(0.5-(a)))

//**************************************************************
//
//	Class name: Point
//
//      9/26/99:  introducing coord array, before Point was a 3D
//      point with the third coordinate 0 if the calculations were
//      done in 2D. The problem now is that the default constructor
//      will need to be in some specified dimension, and this is not
//      good.
//
//      10/22/99: rewriting Point as a template with parameter dim
//      alleviates the above problem and makes it easire to type
//      modifiactions.
//
//      10/24/99: Point<coord_t, dim> suppose to represent points
//      in GEM_Engin. The problem is that to declare
//      Point<coord_t, dim>, dim has to be set, but it cannot be set
//      in GEM_Engin class data members declaration. Of course, I can
//      make dim to be set to the default value 2 or 3 and then
//      use the declaration Point<>, but then I will not be able to
//      change the dimension of the dispaly. So I decided to come back
//      to Point<coord_t> type template, with Point() being a 2D point,
//      but then the assignment operator can overrite its dimension
//      and coordinates to whatever new values we desire. 
//
//      Author: Pawel Gajer
//
//**************************************************************
template <class T = coord_t>
class Point
{
    friend class MesaPlot;
    friend class DrawGraph;
    
    friend ostream &operator<< <T>( ostream &output,
                                    const Point<T> &p );

    //multiplication by a number (made friend for commutativity)
    friend const Point<T>
    operator* <T>( const Point<T> &p, T constant );
    
    friend const Point<T>
    operator* <T>( T constant, const Point<T> &p);

    friend const Point<T>
    operator* <T>( float constant, const Point<T> &p);
    
    friend const Point<T>
    operator/ <T>( const Point<T> &p, T constant );
    
public:
    //default costructor sets dim to 3 and all coordinates to 0
    Point();
    Point(T x, T y); 
    Point(T *_coord, int _dim);
    Point(T x, T y, T z);
    
    const Point<T> &operator=(const Point<T> & rhs); // assign Point
    Point( const Point<T> & init );                  // copy constructor
    ~Point(){ delete [] coord; }                     // destructor
    

    // overloaded arithmetic operators
    const Point<T> operator+(const Point<T> &p) const;
    Point<T> operator+=(const Point<T> &p);
    const Point<T> operator-(const Point<T> &p) const;
    Point<T> operator-() const;
    Point<T> operator-=(const Point<T> &p);
    const T operator*(const Point<T> &p);            // scalar product
    const T operator*(const Point<T> &p) const;      // scalar product
    
    const Point<T> operator*=(T constant);
    const Point<T> operator*=(float constant);
    const Point<T> operator/=(T constant);           // pt <- pt/const
    const Point<T> operator/=(float constant);       // pt <- pt/const
    
    bool operator<(const Point<T> &p) const;
    bool operator==( const Point<T> &right ) const;

    const unsigned long norm() const ; // the norm of the vector associated
                                // with the point
    const float fnorm() const ; // the norm
    const float fnorm2() const ; // the norm 
    const unsigned long norm2() const;      // the square of the norm

    // scalar porduct of v-u and w-u, without the need to compute the differences
    // v-u, w-u first.
    const T scal_pr(const Point<T> &u,
                    const Point<T> &v,
                    const Point<T> &w);
    
   // access operators
    T getX() const { return coord[0]; }
    T getY() const { return coord[1]; }
    T getZ() const { return coord[2]; }

    // modification operators
    void setX(T x){ coord[0] = x; }
    void setY(T y){ coord[1] = y; }
    void setZ(T z){ coord[2] = z; }
    void set_to_zero(){coord[0] = 0;
                       coord[1] = 0;
                       if(dim == 3)
                       coord[2] = 0;
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
//      Point constructors
//
//**************************************************************
template<class T>
Point<T>::Point() :  dim(3), coord(new T[dim])
{
#ifdef DEBUG_POINT
    cout << "Enterting constructor for point of dim " << dim << endl;
#endif
    
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    
#ifdef DEBUG_POINT
    cout << "Leaving constructor for point of dim " << dim << endl;
#endif 
}

template<class T>
Point<T>::Point(T x, T y) :  dim(2), coord(new T[dim])
{
#ifdef DEBUG_POINT
    cout << "Enterting constructor for point of dim " << dim << endl;
#endif
    
    coord[0] = x;
    coord[1] = y;
    
#ifdef DEBUG_POINT
    cout << "Leaving constructor for point of dim " << dim << endl;
#endif 
}

template<class T>
Point<T>::Point(T *_coord, int _dim) : dim(_dim), coord(new T[_dim])
{
#ifdef DEBUG_POINT
    cout << "Enterting 2nd constructor for point of dim " << dim << endl;
#endif
    
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = _coord[ind];
    
#ifdef DEBUG_POINT
    cout << "Leaving 2nd constructor for point of dim " << dim << endl;
#endif 
}

template<class T>
Point<T>::Point(T x, T y, T z) : dim(3), coord(new T[dim])
{
#ifdef DEBUG_POINT
    cout << "Enterting 3nd constructor for point of dim " << dim << endl;
#endif

    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
    
#ifdef DEBUG_POINT
    cout << "Leaving 3nd constructor for point of dim " << dim << endl;
#endif 
}   
//**************************************************************
//
//	Method name : Copy constructor
//
//**************************************************************
template<class T>
Point<T>::Point( const Point<T> &init ) :
        dim(init.dim), coord(new T[dim])
{
#ifdef DEBUG_POINT
    cout << "enterting copy constructor for point of dim "
         << dim << endl;
#endif
   for ( size_t i = 0; i < dim; i++ )
       coord[ i ] = init.coord[ i ];  // copy init into object

#ifdef DEBUG_POINT
    cout << "Leaving copy constructor for Point of dim "
         << dim << endl;
#endif
}

//**************************************************************
//
//	Method name : assignment operator (=)
//
//**************************************************************
template<class T>
const Point<T> &Point<T>::operator=( const Point<T> &rhs )
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
//	Description : Output Point with overloaded stream
//      insertion operator.
//
//**************************************************************
template<class T>
ostream &operator<<( ostream &output, const Point<T> &p )
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
const Point<T> Point<T>::operator+(const Point<T> &p) const
{
//#define DEBUG_POINT
#ifdef DEBUG_POINT
    cout << "Enterting (+) for points  " << (*this)
         << " and " << p << endl;
#endif  
    T *new_coord = new T[dim];
    for(size_t ind = 0; ind < dim; ind++)
        new_coord[ind] = coord[ind] + p.coord[ind];

    Point<T> result(new_coord, dim);
    delete [] new_coord;

#ifdef DEBUG_POINT
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
Point<T> Point<T>::operator+=(const Point<T> &p)
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
const Point<T> Point<T>::operator-(const Point<T> &p) const
{
    T *new_coord = new T[dim];
    for(size_t ind = 0; ind < dim; ind++)
        new_coord[ind] = coord[ind] - p.coord[ind];

    Point result(new_coord, dim);
    delete [] new_coord;

    return result;
}
// negative
template<class T>
Point<T> Point<T>::operator-() const
{
    return (*this) * (-1);
}

//**************************************************************
//
//	Method name : subtraction assignment (-=)
//
//**************************************************************
template<class T>
Point<T> Point<T>::operator-=(const Point<T> &p)
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
const Point<T> operator*( T constant, const Point<T> &p)
{
    T *new_coord = new T[p.dim];
    
    for(size_t ind = 0; ind < p.dim; ind++)
        new_coord[ind] = constant * p.coord[ind];

    Point<T> result(new_coord, p.dim);
    delete [] new_coord;

    return result;
}

template<class T>
const Point<T> operator*( float constant, const Point<T> &p)
{
    T *new_coord = new T[p.dim];
    
    for(size_t ind = 0; ind < p.dim; ind++)
        new_coord[ind] = (T)(ROUND(constant * p.coord[ind]));

    Point<T> result(new_coord, p.dim);
    delete [] new_coord;

    return result;
}

template<class T>
const Point<T> operator*( const Point<T> &p,
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
const Point<T> Point<T>::operator*=( T constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] *= constant;
    
    return *this;
}

template<class T>
const Point<T> Point<T>::operator*=( float constant )
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
const Point<T> Point<T>::operator/=( T constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] /= constant;
    
    return *this;
}

template<class T>
const Point<T> Point<T>::operator/=( float constant )
{
    for(size_t ind = 0; ind < dim; ind++)
        coord[ind] = (T)((float)coord[ind]/constant);
    
    return *this;
}

//**************************************************************
//
//	Method name : (/) 
//
//      Description :  pt/constant
//
//**************************************************************
template<class T>
const Point<T> operator/( const Point<T> &p, T constant )
{
    T *new_coord = new T[p.dim];
    
    for(size_t ind = 0; ind < p.dim; ind++)
        new_coord[ind] = p.coord[ind] / constant;
    
    Point<T> result(new_coord, p.dim);
    delete [] new_coord;

    return result;
}

//**************************************************************
//
//	Method name : (*) scalar product 
//
//      Description : scalar product of vectors associated with
//      points
//
//**************************************************************
template<class T>
const T Point<T>::operator*(const Point<T> &p)
{
    T result = 0;
    for(size_t ind = 0; ind < p.dim; ind++)
        result += coord[ind] * p.coord[ind];

    return result;
}

template<class T>
const T Point<T>::operator*(const Point<T> &p) const
{
    T result = 0;
    for(size_t ind = 0; ind < p.dim; ind++)
        result += coord[ind] * p.coord[ind];

    return result;
}

//**************************************************************
//
//       scal_pr()
//       scalar porduct of uv = v-u and uw = w-u,
//       without the need to compute the differences v-u, w-u.
//
//**************************************************************
template<class T>
const T scal_pr(const Point<T> &u,
                const Point<T> &v,
                const Point<T> &w)
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
bool Point<T>::operator<(const Point<T> &p) const
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
const unsigned long Point<T>::norm() const
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
const float Point<T>::fnorm() const
{
    return sqrt( norm2() );
}

//**************************************************************
//
//	Method name : norm2
//
//	Description : the square of the norm
//
//**************************************************************
template<class T>
const unsigned long Point<T>::norm2() const
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
const float Point<T>::fnorm2() const
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
bool Point<T>::operator==( const Point<T> &right ) const
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
