#ifndef UBLAS_EXT_CONFIG_H
#define UBLAS_EXT_CONFIG_H

#define BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
#define BOOST_UBLAS_EXT_ENABLE_MACROS

#ifdef BOOST_UBLAS_EXT_ENABLE_MACROS
#	define END last_elem<>()
#	define ALL all_elem()
#	define REV_ALL rev_all()
#	define RG(x, y)		(boost::numeric::ublas::basic_slice_ext<>((x),(y)))
#	define SL(x, y, z)  (boost::numeric::ublas::basic_slice_ext<>((x),(y),(z)))
#endif

#endif //UBLAS_EXT_CONFIG_H
