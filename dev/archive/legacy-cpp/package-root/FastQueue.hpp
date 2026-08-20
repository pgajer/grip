// FastQueue.hpp a header file for FastQueue class

#ifndef FAST_QUEUE_HPP
#define FAST_QUEUE_HPP

#include <iostream>
#include <cassert>

#include "Debug.hpp"

typedef unsigned short size_tt;


//**************************************************************
//
//	Class name : FastQueue
//
//	Description : queue data structure based on an array
//      not a circular array. We need to know how many elements
//      will be entered into the queue.
//
//      Date: May 9, 2000
//
//**************************************************************
template <class Item>
class FastQueue
{
    template <class T>
    friend std::ostream& operator<<(std::ostream&, const FastQueue<T>&);

public:
    // CONSTRUCTORS, ASSIGNMENT OP, and DESTRUCTOR
    FastQueue( unsigned long _size = 1 );
    FastQueue(const FastQueue& source);           // copy constructor
    FastQueue &operator =(const FastQueue& source); // assignment operator
    ~FastQueue();
    
    // MODIFICATION functions
    void enqueue(const Item& entry);
    Item dequeue(); // return the first element and remove it from the queue
    Item top(){ return ptr[frontPtr]; } // return the first element
    void empty(){ frontPtr = 0; backPtr = 0; count = 0;}
    // CONSTANT functions
    unsigned long get_count( ) const { return count; }
    unsigned long get_size( ) const { return size; }
    bool is_empty( ) const { return (count == 0); }

private:
    Item *ptr;        // an array holding queue elements
    unsigned long size;     // size of the array
    unsigned long count;    // Total number of items in the queue
    unsigned long frontPtr; // indices pointing to the front and back
    unsigned long backPtr;  // of the queue
};

//==============================================================
//
//          TEMPLATE CLASS MEMBER FUNCTION DEFINITIONS
//
//==============================================================

//**************************************************************
//
//	class constructor 
//
//**************************************************************
template <class Item>
FastQueue<Item>::FastQueue(unsigned long _size)
        : ptr(new Item[_size]),
          size(_size),
          count(0),
          frontPtr(0),
          backPtr(0)
{
}

//**************************************************************
//
//     Copy constructor
//
//**************************************************************
template< class Item >
FastQueue< Item >::FastQueue( const FastQueue< Item > &init )
        : size( init.size )
{
#if DEBUG_FQ
    debug("coping an object of size " << size);
    for (unsigned long i = 0; i < size; i++ )
        cout << ptr[ i ] << ' ';
    cout << endl;
#endif
   ptr = new Item[ size ]; // create space for FastQueue
   assert( ptr != 0 );    // terminate if memory not allocated

   for (unsigned long i = 0; i < init.size; i++ )
       ptr[ i ] = init.ptr[ i ];  // copy init into object
}

//**************************************************************
//
//	class destructor
//
//**************************************************************
template< class T >
FastQueue< T >::~FastQueue()
{
#define DEBUG_FQ 0
#if DEBUG_FQ 
    debug("destructing an object of size " << size);
    for (unsigned long i = 0; i < size; i++ )
        cout << ptr[ i ] << ' ';
    cout << endl;
#endif
//    debug("before delete [] ptr in FastQueue; size="<<size<<", count="<<count);
    
    delete [] ptr;            // reclaim space for FastQueue

//    debug("FastQueue destroyed");
        
}

//**************************************************************
//
//	Overloaded assignment operator (=)
//
//**************************************************************
template< class Item >
FastQueue< Item > &FastQueue< Item >::operator=( const FastQueue< Item > &right )
        // const return avoids: ( a1 = a2 ) = a
{
    if ( &right != this ) // check for self-assignment
    {  
        // for vectors of different sizes, deallocate original
        // left side vector, then allocate new left side vector.
        if ( size != right.size )
        {
            delete [] ptr;          // reclaim space
            size = right.size;      // resize this object
            ptr = new Item[ size ];// create space for vector copy
        }
        
        for ( unsigned long i = 0; i < right.size; i++ )
            ptr[ i ] = right.ptr[ i ];  // copy FastQueue into object
    }
    
    return *this;   // enables x = y = z;
}

//**************************************************************
//
//	Overloaded output operator (<<)
//
//**************************************************************
template< class Item >
std::ostream &operator<<(std::ostream &output,
                         const FastQueue< Item > &a )
{
    if(a.count){
        for (unsigned long i = a.frontPtr; i < a.backPtr; i++ )
            output << a.ptr[ i ] << ' ';
        output << std::endl;
    } else
        output << "FastQueue is empty" << std::endl;
    
    return output;   // enables cout << x << y;
}

//**************************************************************
//
//	enqueue()
//
//**************************************************************
template< class Item >
void FastQueue< Item >::enqueue(const Item &entry)
{
    assert(count < size);
    count++;
    ptr[backPtr++] = entry;
}

//**************************************************************
//
//	dequeue()
//
//**************************************************************
template< class Item >
Item FastQueue< Item >::dequeue()
{
    assert(count > 0 && frontPtr+1 < size);
    count--;
    return ptr[frontPtr++];
}

#endif
