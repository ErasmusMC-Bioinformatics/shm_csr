/* script adapted from https://www.nu42.com/2021/07/windows-c-time-in-nanoseconds.html */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double
time_now(void)
{
    struct timespec ts;

    if (timespec_get(&ts, TIME_UTC) != TIME_UTC)
    {
        fputs("timespec_get failed!", stderr);
        return 0;
    }
    return (double)ts.tv_sec + ((double)ts.tv_nsec / (double)1000000000);
}

int main(void)
{
    printf("%lf\n", time_now());
    return EXIT_SUCCESS;
}
