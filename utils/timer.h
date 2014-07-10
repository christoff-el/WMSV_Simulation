//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
// Simulation of Multidimensional Stochastic //
//            Volatility Models              //
//                                           //
//           Implementation Code             //
//                                           //
//         Christopher Davis - 2014          //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/*

Description: Simple timer.
             
*/

#ifndef TIMER_H
#define TIMER_H 1

#include <iostream>
#include <string>
#include <sys/time.h>

class Timer{

private:
    timeval _start;
    timeval _stop;
    
public:
    
    //Start the timer:
    void
    start(){
        gettimeofday(&_start, 0);
    }
    
    //Stop the timer:
    void
    stop(){
        gettimeofday(&_stop, 0);
    }
    
    //Get the elapsed time:
    double
    elapsed(){
        long long int usecs_start = _start.tv_sec * 1000000. + _start.tv_usec;
        long long int usecs_stop =  _stop.tv_sec * 1000000. + _stop.tv_usec;
        return double(usecs_stop - usecs_start) / 1000000.;
    }
    		
};


#endif // TIMER_H