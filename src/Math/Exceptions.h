/*! \file Exceptions.h
 *  \brief header file that defines exceptions
 */
// Function.h

#ifndef NBODYEXCEPTION_H
#define NBODYEXCEPTION_H

#include <exception>
#include <string>
#include <gsl/gsl_errno.h>

namespace Math
{

class gsl_error : public std::exception {

private:
	int gsl_errno;
	std::string errmsg;

public:
	gsl_error(int gsl_errno) :
		gsl_errno(gsl_errno),
		errmsg(gsl_strerror(gsl_errno))
	{
	}

	int get_gsl_errno() {
		return gsl_errno;
	}

	const char *what() const noexcept
	{
		return errmsg.c_str();
	}

};

}


#endif
