/* script adapted from https://www.nu42.com/2021/07/windows-c-time-in-nanoseconds.html */
#include <stdio.h>
#include <time.h>

int main(void)
{
    struct timespec ts;

    if (timespec_get(&ts, TIME_UTC) != TIME_UTC)
    {
        fputs("timespec_get failed!", stderr);
        return 1;
    }
    printf("%d.%d\n", ts.tv_sec, ts.tv_nsec);
    return 0;
}
