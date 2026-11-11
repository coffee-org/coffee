/**
 * @file    coffee_compat.h
 * @brief   Backward-compatibility shims for coffee optical modules
 *
 * Provides thin inline shims for removed semaphore creation calls and
 * OpticsMaterials naming aliases when building against modern milk.
 */

#ifndef COFFEE_COMPAT_H
#define COFFEE_COMPAT_H

#include <string.h>

#include "milk_runnable.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "OpticsMaterials/OpticsMaterials.h"

/**
 * COREMOD_MEMORY_image_set_createsem - Legacy no-op semaphore initializer
 * @IDname: Image stream name (unused)
 * @NBsem:  Requested semaphore count (unused)
 *
 * Semaphores are created automatically by ImageStreamIO on stream creation.
 *
 * Return: 0
 */
static inline imageID COREMOD_MEMORY_image_set_createsem(
    const char *IDname __attribute__((unused)),
    long        NBsem __attribute__((unused)))
{
    return 0;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846264338328
#endif

#define COFFEE_MATCODE_CSTPHA 98

/**
 * OpticsMaterials_n - Compute refractive index at wavelength
 * @material: Material identifier code
 * @lambda:   Wavelength in meters
 *
 * Return: Refractive index n(lambda).
 */
static inline double OpticsMaterials_n(
    int    material,
    double lambda)
{
    if (material == COFFEE_MATCODE_CSTPHA)
    {
        return 1.0 + lambda * 1.0e6;
    }
    return OPTICSMATERIALS_n(material, lambda);
}

/**
 * OpticsMaterials_code - Resolve material identifier code from string name
 * @name: Material name string
 *
 * Return: Numeric material code, or -1 if unknown.
 */
static inline int OpticsMaterials_code(
    const char *name)
{
    if (name != NULL && strcmp(name, "CstPha") == 0)
    {
        return COFFEE_MATCODE_CSTPHA;
    }
    return OPTICSMATERIALS_code(name);
}

/**
 * OpticsMaterials_name - Resolve material string name from identifier code
 * @code: Numeric material code
 *
 * Return: Material name string, or NULL if unknown.
 */
static inline const char *OpticsMaterials_name(
    int code)
{
    if (code == COFFEE_MATCODE_CSTPHA)
    {
        return "CstPha";
    }
    return OPTICSMATERIALS_name(code);
}

/**
 * OpticsMaterials_pha_lambda - Compute optical phase shift across thickness
 * @material: Material identifier code
 * @z:        Physical thickness in meters
 * @lambda:   Wavelength in meters
 *
 * Return: Phase shift in radians.
 */
static inline double OpticsMaterials_pha_lambda(
    int    material,
    double z,
    double lambda)
{
    if (material == COFFEE_MATCODE_CSTPHA)
    {
        return 2.0 * M_PI * z * 1.0e6;
    }
    return OPTICSMATERIALS_pha_lambda(material, z, lambda);
}

/**
 * coffee_dcvar_f - Read numeric CLI variable as double regardless of union tag
 * @id: Variable table index from variable_ID()
 *
 * In milk's @calc evaluator, integer literals create type 1 (long, value.l)
 * while floating-point literals create type 0 (double, value.f). This helper
 * safely reads either numeric representation as double.
 *
 * Return: Numeric value of variable, or 0.0 if @id < 0 or non-numeric.
 */
static inline double coffee_dcvar_f(
    variableID id)
{
    if (id < 0)
    {
        return 0.0;
    }
    if (dcvar[id].type == 1)
    {
        return (double) dcvar[id].value.l;
    }
    if (dcvar[id].type == 0)
    {
        return dcvar[id].value.f;
    }
    return 0.0;
}

#endif /* COFFEE_COMPAT_H */
