#ifndef __DEFINES_H__
#define __DEFINES_H__

// Period measured in milliseconds
#define MIN_TI 10
#define MAX_TI 1000
#define MEAN_TI 500
#define STDEV_TI 200

// Parameters for short critical section (microseconds)
// Mean and standard deviation used only with normal distributions
#define MIN_SHORT_CSLEN 1
#define MAX_SHORT_CSLEN 15
#define MEAN_SHORT_CSLEN 8
#define STDEV_SHORT_CSLEN 4

// Parameters for moderate critical section (microseconds)
#define MIN_MODERATE_CSLEN 1
#define MAX_MODERATE_CSLEN 100
#define MEAN_LONG_CSLEN 50
#define STDEV_LONG_CSLEN 30

// Parameters for long critical section (microseconds)
#define MIN_LONG_CSLEN 5
#define MAX_LONG_CSLEN 1280

// Type of critical section (short or moderate or long)
typedef enum {
	SHORT = 0, 
	MODERATE, 
	LONG,
} CriticalDuration;

// Type of spin locks
typedef enum {
	FIFO = 0, 
	PRIO_UNORDERED, 
	PRIO_FIFO,
} SpinlockType;

typedef unsigned int ResourceID;
typedef unsigned int TaskID;

#endif // __DEFINES_H__
