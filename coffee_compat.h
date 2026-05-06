/**
 * @file    coffee_compat.h
 * @brief   Backward-compatibility shims for legacy
 *          coffee code
 *
 * Provides stubs for removed functions and extern
 * declarations for functions missing headers.
 *
 * Include this header AFTER CLIcore.h in every
 * coffee .c file that uses the old API.
 */

#ifndef COFFEE_COMPAT_H
#define COFFEE_COMPAT_H

#include "CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"


/* ============================================
 * COREMOD_MEMORY_image_set_createsem:
 * Removed — semaphores are now created
 * automatically by ImageStreamIO.
 * Stub to no-op.
 * ============================================ */

static inline imageID
COREMOD_MEMORY_image_set_createsem(
    const char *IDname __attribute__((unused)),
    long        NBsem  __attribute__((unused)))
{
    return 0;
}


/* ============================================
 * waitforsemID:
 * Exists in stream_sem.c but has no header.
 * Provide declaration.
 * ============================================ */

extern void *waitforsemID(void *ID);


/* ============================================
 * COREMOD_MEMORY_image_set_status:
 * May have different signature — provide
 * wrapper if needed.
 * ============================================ */


#endif /* COFFEE_COMPAT_H */
