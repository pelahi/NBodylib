/*! \file Exceptions.cxx
 *  \brief subroutines for handling exceptions
 */

#include <Exceptions.h>

using namespace std;
namespace Math
{

void throw_exception_gsl_handler(const char *reason, const char *file, int line, int gsl_errno)
{
    throw Math::gsl_error(reason, file, line, gsl_errno, gsl_strerror(gsl_errno));
}
void install_gsl_error_handler()
{
    gsl_set_error_handler(&throw_exception_gsl_handler);
}


}
