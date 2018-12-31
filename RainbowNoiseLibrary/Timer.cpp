/*
 * Timer.cpp
 *
 * Li-Yi Wei
 * 4/15/2003
 *
 */

#include "Timer.hpp"

Timer::Timer(void) : _started(0)
{
    // nothing to do
}

Timer::~Timer(void)
{
    // nothing to do
}

int Timer::Start(void)
{
#ifdef WIN32
    if(timeBeginPeriod(1) == TIMERR_NOERROR)
    {
        _started = 1;
        _startRealTime = timeGetTime();
        return 1;
    }
#else
    struct timeval t;
    if (gettimeofday(&t, NULL) == 0)
    {
        _started = 1;
        _startRealTime = (t.tv_sec + t.tv_usec / 1000000.0);
        return 1;
    }
#endif
    else
    {
        return 0;
    }
}

int Timer::Stop(void)
{
    _started = 0;
#ifdef WIN32
    _stopRealTime = timeGetTime();
    timeEndPeriod(1);
#else
    struct timeval t;
    if (gettimeofday(&t, NULL) == 0)
    {
        _stopRealTime = (t.tv_sec + t.tv_usec / 1000000.0);
        return 1;
    }
    else
    {
        return 0;
    }
#endif
    return 1;
}

double Timer::ElapsedTime(const Time & startTime,
			  const Time & endTime) const
{
#ifdef WIN32
    return (_stopRealTime - _startRealTime)*0.001;
#else
    return (_stopRealTime - _startRealTime);
#endif
}

double Timer::ElapsedTime(void) const
{
    return ElapsedTime(_startRealTime, _stopRealTime);
}

double Timer::CurrentTime(void) const
{
    if(_started)
    {
        return 0;
    }
    else
    {
#ifdef WIN32
        if(timeBeginPeriod(1) == TIMERR_NOERROR)
        {
            Time value = timeGetTime();
            timeEndPeriod(1);
            return value * 0.001;
        }
#else
        struct timeval t;
        if (gettimeofday(&t, NULL) == 0)
        {
            return (t.tv_sec + t.tv_usec / 1000000.0);
        }
#endif
        else
        {
            return 0;
        }
    }
}
