// Fixed-capacity array-backed queue.

#ifndef FAST_QUEUE_HPP
#define FAST_QUEUE_HPP

#include <cstdint>
#include <stdexcept>

using size_tt = uint32_t;


//**************************************************************
//
//	Class name : FastQueue
//
//	Description : queue data structure based on an array
//      not a circular array. We need to know how many elements
//      will be entered into the queue.
//
//**************************************************************
template <class Item>
class FastQueue
{
public:
    // CONSTRUCTORS, ASSIGNMENT OP, and DESTRUCTOR
    FastQueue( unsigned long _size = 1 );
    FastQueue(const FastQueue& source);           // copy constructor
    FastQueue &operator =(const FastQueue& source); // assignment operator
    ~FastQueue();
    
    // MODIFICATION functions
    void enqueue(const Item& entry);
    Item dequeue(); // return the first element and remove it from the queue
    Item top(){
        if (count == 0 || frontPtr >= backPtr || frontPtr >= size)
            throw std::underflow_error("cannot read from an empty FastQueue");
        return ptr[frontPtr];
    } // return the first element
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
        : ptr(new Item[init.size]),
          size(init.size),
          count(init.count),
          frontPtr(init.frontPtr),
          backPtr(init.backPtr)
{
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
    
    delete [] ptr;            // reclaim space for FastQueue

        
}

//**************************************************************
//
//	Overloaded assignment operator (=)
//
//**************************************************************
template< class Item >
FastQueue< Item > &FastQueue< Item >::operator=( const FastQueue< Item > &right )
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
        count = right.count;
        frontPtr = right.frontPtr;
        backPtr = right.backPtr;
    }
    
    return *this;   // enables x = y = z;
}

//**************************************************************
//
//	enqueue()
//
//**************************************************************
template< class Item >
void FastQueue< Item >::enqueue(const Item &entry)
{
    if (count >= size || backPtr >= size)
        throw std::overflow_error("FastQueue capacity exceeded");
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
    if (count == 0 || frontPtr >= backPtr || frontPtr >= size)
        throw std::underflow_error("cannot dequeue from an empty FastQueue");
    count--;
    return ptr[frontPtr++];
}

#endif
