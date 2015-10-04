#ifndef TIME__H
#define TIME__H

#define CLOCK_MONOTONIC 0

#include <sys/time.h>
#include <time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>


typedef int clockid_t;

static int
clock_gettime(clockid_t clk_id, struct timespec *tp)
{
	struct timespec ts;
	clock_serv_t cclock;
	mach_timespec_t mts;
	host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
	clock_get_time(cclock, &mts);
	mach_port_deallocate(mach_task_self(), cclock);
	ts.tv_sec = mts.tv_sec;
	ts.tv_nsec = mts.tv_nsec;
	*tp = ts;
	return 0;
}
#endif

#endif
