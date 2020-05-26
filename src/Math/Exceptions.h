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
	std::string reason;
	std::string file;
	int line;
	int gsl_errno;
	std::string errmsg;

public:
	gsl_error(const char *reason, const char *file, int line, int gsl_errno, const char *errmsg) :
		reason(reason),
		file(file),
		line(line),
		gsl_errno(gsl_errno),
		errmsg(errmsg)
	{
	}

	std::string get_reason() {
		return reason;
	}

	std::string get_file() {
		return file;
	}

	int get_line() {
		return line;
	}

	int get_gsl_errno() {
		return gsl_errno;
	}

    std::string get_errmsg() {
        return errmsg;
    }

};

void throw_exception_gsl_handler(const char *reason, const char *file, int line, int gsl_errno);
void install_gsl_error_handler(); 
}


#endif
