/* declared in src/include/Defn.h 
   and
   defined in src/main/attrib.c */

#include <Rversion.h>
#if R_VERSION < R_Version(4, 6, 0)
# define R_class(x) R_data_class(x, FALSE)
#endif

SEXP R_data_class(SEXP obj, Rboolean singleString);
