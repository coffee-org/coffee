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
 * OpticsMaterials shims:
 * Maps legacy camelCase OpticsMaterials_* to
 * current OPTICSMATERIALS_* API.
 * ============================================ */

#include "OpticsMaterials/OpticsMaterials.h"

static inline double OpticsMaterials_n(int material, double lambda)
{
    return OPTICSMATERIALS_n(material, lambda);
}

static inline int OpticsMaterials_code(const char *name)
{
    return OPTICSMATERIALS_code(name);
}

static inline const char *OpticsMaterials_name(int code)
{
    return OPTICSMATERIALS_name(code);
}

static inline double OpticsMaterials_pha_lambda(int material, double z, double lambda)
{
    return OPTICSMATERIALS_pha_lambda(material, z, lambda);
}


/* ============================================
 * INSERT_STD_CLIfunction:
 * Compatibility macro for legacy CLI function
 * argument array parsing and execution.
 * ============================================ */

extern errno_t CLI_checkarg_array(CLICMDARGDEF fpscliarg[], int nbarg);
extern void *get_farg_ptr(char *tag, long *fpsi);

#ifndef STD_FARG_LINKfunction
#define STD_FARG_LINKfunction                                                  \
    for (int argi = 0; argi < (int) (sizeof(farg) / sizeof(CLICMDARGDEF));     \
         argi++)                                                               \
    {                                                                          \
        long  fpsi           = -1;                                             \
        void *ptr            = get_farg_ptr(farg[argi].fpstag, &fpsi);         \
        *(farg[argi].valptr) = ptr;                                            \
        if (farg[argi].indexptr != NULL)                                       \
        {                                                                      \
            *(farg[argi].indexptr) = fpsi;                                     \
        }                                                                      \
    }
#endif

#ifndef INSERT_STD_CLIfunction
#define INSERT_STD_CLIfunction                                                 \
    static errno_t CLIfunction(void)                                           \
    {                                                                          \
        errno_t retval = CLI_checkarg_array(farg, CLIcmddata.nbarg);           \
        if (retval == RETURN_SUCCESS)                                          \
        {                                                                      \
            STD_FARG_LINKfunction return compute_function();                   \
        }                                                                      \
        if (retval == RETURN_CLICHECKARGARRAY_HELP)                            \
        {                                                                      \
            return RETURN_SUCCESS;                                             \
        }                                                                      \
        if (retval == RETURN_CLICHECKARGARRAY_FUNCPARAMSET)                    \
        {                                                                      \
            return RETURN_SUCCESS;                                             \
        }                                                                      \
        return retval;                                                         \
    }
#endif

#endif /* COFFEE_COMPAT_H */
