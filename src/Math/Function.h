/*! \file Function.h
 *  \brief header file that defines several function structures
 */
// Function.h

#ifndef NBODYFUNCTION_H
#define NBODYFUNCTION_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace Math
{

/// \struct Math::gsl_math_function
/// \brief gsl function wrapper
struct gsl_math_function
{
  int (* function) (const gsl_vector *params, void * data, gsl_vector *output);
};

/// \struct Math::gsl_math_function
/// \brief gsl function wrapper differentiating wrt params, storing info in Jacobian matrix
struct gsl_math_function_df
{
  int (* function) (const gsl_vector * params, void * data, gsl_matrix * outputJacobian);
};

/// \struct Math::gsl_fitting_data
struct gsl_fitting_data {
  size_t n;
  double * x;
  double * y;
  size_t npar;
  int * iparindex;
  size_t nallparams;
  double *allparams;
};

/// \struct Math::math_function
/// \brief single dependent variable function
struct math_function
{
  Double_t (* function) (Double_t x, void * params);
  int (* gsl_function) (const gsl_vector *params, void * data, gsl_vector *output);
  int (* gsl_function_df) (const gsl_vector *params, void * data, gsl_matrix * outputJacobian);
  void * params;
};

/// \struct Math::math_multidim_function
/// \brief multi-dimensional function
struct math_multidim_function
{
  Double_t (* function) (Double_t *x, int ndim, void * params);
  int ndim;
  void * params;
};
}
#endif
