/**
 * @file debug.h
 * @brief Defines checks and asserts for program testing and debugging.
 */

#ifndef RXI_DEBUG_H
#define RXI_DEBUG_H

#include <assert.h>

/// @brief Check condition and abort on fail.
#define ASSERT(condition) assert(condition)

/// @brief Check condition and print message on fail without abort.
#ifndef NDEBUG
# define CHECK(condition)                                                     \
    do                                                                        \
      {                                                                       \
        if (!condition)                                                       \
          {                                                                   \
            printf ("%s:%u: %s: ", __FILE__, __LINE__, __PRETTY_FUNCTION__);  \
            printf ("Check `%s' failed.\n", #condition);                      \
          }                                                                   \
      }                                                                       \
    while(0)
#else
# define CHECK(condition) ((void)0)
#endif

/// @brief Using this instead of `printf()` to print debug messages in stdio.
#ifndef NDEBUG
# define DEBUG(...)                                                           \
    do                                                                        \
      {                                                                       \
        printf ("%s | ", __TIME__);                                           \
        printf ("%s:%u: %s: ", __FILE__, __LINE__, __func__);                 \
        printf (__VA_ARGS__);                                                 \
        printf ("\n");                                                        \
      }                                                                       \
    while(0)
#else
# define DEBUG(...) ((void)0)
#endif

#endif // RXI_DEBUG_H
