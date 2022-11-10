#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void kalmanLoopC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kalmanLoopTempC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mleLoopCoe(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mleLoopCpe(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"kalmanLoopC",     (DL_FUNC) &kalmanLoopC,     14},
    {"kalmanLoopTempC", (DL_FUNC) &kalmanLoopTempC, 10},
    {"mleLoopCoe",      (DL_FUNC) &mleLoopCoe,      10},
    {"mleLoopCpe",      (DL_FUNC) &mleLoopCpe,      10},
    {NULL, NULL, 0}
};

void R_init_LakeMetabolizer(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
