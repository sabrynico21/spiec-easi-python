### configure.ac					-*- Autoconf -*-

AC_PREREQ(2.62)

AC_INIT([huge],2.0.1)
AC_PATH_PROG([RBIN], [R])

OPENMP_CXXFLAGS=""

if test `uname -s` = "Linux"
then
  OPENMP_CXXFLAGS="\$(SHLIB_OPENMP_CXXFLAGS)"
fi

if test `uname -s` = "Darwin"
then
  OPENMP_CXXFLAGS="\$(SHLIB_OPENMP_CXXFLAGS)"
  ac_pkg_openmp=no
  AC_MSG_CHECKING([whether OpenMP will work in a package])
  AC_LANG_PUSH([C++])
  AC_LANG_CONFTEST(
  [AC_LANG_SOURCE([[
  #include <omp.h>
  #ifndef _OPENMP
  #error "OpenMP compilation is not enabled"
  #endif
  extern "C" void huge_openmp_probe(int *result) {
    int total = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:total)
    for(int i = 0; i < 4; ++i)
      total += i;
    result[0] = total;
    result[1] = omp_get_max_threads();
  }
  ]])])
  PKG_CXXFLAGS="${OPENMP_CXXFLAGS}" PKG_LIBS="${OPENMP_CXXFLAGS}" "$RBIN" CMD SHLIB conftest.cpp 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "$RBIN" --vanilla -q -e "so <- paste0('conftest', .Platform\$dynlib.ext); dyn.load(so); probe <- unname(unlist(.C('huge_openmp_probe', integer(2)))); stopifnot(probe[[1]] == 6L, probe[[2]] >= 1L); dyn.unload(so)" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
  AC_LANG_POP([C++])
  AC_MSG_RESULT([${ac_pkg_openmp}])
  if test "${ac_pkg_openmp}" = no; then
    OPENMP_CXXFLAGS=''
  fi
fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
