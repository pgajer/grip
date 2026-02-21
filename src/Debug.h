#ifndef DEBUG_H_
#define DEBUG_H_

#include <iostream>

#ifdef GRIP_DEBUG
#define debug(a) std::cout << "<"<< __FILE__ << ", " << __FUNCTION__ << ">[" << __LINE__ << "] " << a << std::endl;
#else
#define debug(a) do {} while(0)
#endif

#endif // DEBUG_H_
