#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([LAMP], [3.0.0], [lamp_staff@googlegroups.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([LampCore.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_RANLIB
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

# Checks for libraries.
GTEST_LIB_CHECK([1.7.0], [true], [true])

BOOST_REQUIRE([1.41.0])
BOOST_FILESYSTEM
BOOST_IOSTREAMS
BOOST_PROGRAM_OPTIONS
BOOST_SYSTEM

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([dup2 pow sqrt strchr])

AC_CONFIG_FILES([Makefile
                 frepattern/Makefile
                 functions/Makefile
                 lcm53/Makefile
                 flowers/Makefile])
AC_OUTPUT

