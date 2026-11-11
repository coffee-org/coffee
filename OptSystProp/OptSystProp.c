/**
 * @file    OptSystProp.c
 * @brief   Optical system Fresnel propagation engine module descriptor
 */

#define MODULE_SHORTNAME_DEFAULT "optsystprop"
#define MODULE_DESCRIPTION       "Optical System Fresnel Propagation Engine"

#include "milk_runnable.h"
#include "coffee_compat.h"
#include "OptSystProp/OptSystProp.h"

MODULE_DEPS("milkfft", "milkOpticsMaterials", "milkWFpropagate");

/**
 * init_OptSystProp - Initialize OptSystProp module
 *
 * Return: RETURN_SUCCESS
 */
errno_t init_OptSystProp(void)
{
    return RETURN_SUCCESS;
}

static const MILK_RUNNABLE *const runnables[] = {
    NULL,
};

MILK_MODULE_RUNNABLES(coffeeOptSystProp, _module_deps, runnables);
