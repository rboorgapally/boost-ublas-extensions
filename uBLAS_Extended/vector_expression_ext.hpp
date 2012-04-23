#ifndef VECTOR_EXPRESSION_EXT_H
#define VECTOR_EXPRESSION_EXT_H

#include "functional_ext.hpp"
#include <boost/numeric/ublas/vector_expression.hpp>

namespace boost { namespace numeric { namespace ublas {

	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary<E,
		scalar_log_minimum <typename E::value_type> >
		>::result_type
		log_min (const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				scalar_log_minimum <typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		scalar_log_maximum <typename E::value_type> >
		>::result_type
		log_max (const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				scalar_log_maximum <typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E> 
	BOOST_UBLAS_INLINE 
		typename vector_scalar_unary_traits<E, 
		general_vector_scalar_unary< E,
		scalar_minimum<typename E::value_type> > 
		>::result_type 
		min (const vector_expression<E> &e) { 
			typedef typename vector_scalar_unary_traits<E, 
				general_vector_scalar_unary< E, 
				scalar_minimum<typename E::value_type> > 
			>::expression_type expression_type; 
			return expression_type (e ()); 
	} 
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		scalar_maximum<typename E::value_type> >
		>::result_type
		max (const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				scalar_maximum<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}

	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		abs_min_vector<typename E::value_type> >
		>::result_type
		abs_min(const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				abs_min_vector<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		abs_max_vector<typename E::value_type> >
		>::result_type
		abs_max(const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				abs_max_vector<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		mean_vector<typename E::value_type> >
		>::result_type
		mean(const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				mean_vector<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		abs_mean_vector<typename E::value_type> >
		>::result_type
		abs_mean(const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				abs_mean_vector<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_scalar_unary_traits<E,
		general_vector_scalar_unary< E,
		rms_vector<typename E::value_type> >
		>::result_type
		rms(const vector_expression<E> &e) {
			typedef typename vector_scalar_unary_traits<E,
				general_vector_scalar_unary< E,
				rms_vector<typename E::value_type> >
			>::expression_type expression_type;
			return expression_type (e ());
	}
	// (t + v) [i] = t + v [i]
	template<class T1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::result_type
		operator + (const T1 &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (t - v) [i] = t - v [i]
	template<class T1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<const T1, E2, scalar_minus<T1, typename E2::value_type> >::result_type
		operator - (const T1 &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<const T1, E2, scalar_minus<T1, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v + t) [i] = v [i] + t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::result_type
		operator + (const vector_expression<E1> &e1, const T2 &e2) {
			typedef typename vector_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (v - t) [i] = v [i] - t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, const T2, scalar_minus<typename E1::value_type, T2> >::result_type
		operator - (const vector_expression<E1> &e1, const T2 &e2) {
			typedef typename vector_binary_scalar2_traits<E1, const T2, scalar_minus<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
											

	// (v ^ t) [i] = v [i] ^ t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::result_type
		operator ^ (const vector_expression<E1> &e1, const T2 &e2) {
			typedef typename vector_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (abs v) [i] = abs (v [i])
	template<class E> 
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_abs<typename E::value_type> >::result_type
		abs (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_abs<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (exp v) [i] = exp (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_exp<typename E::value_type> >::result_type
		exp (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_exp<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (log v) [i] = log (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_log<typename E::value_type> >::result_type
		log (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_log<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (log10 v) [i] = log10 (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_log10<typename E::value_type> >::result_type
		log10 (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_log10<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (sqrt v) [i] = sqrt (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_sqrt<typename E::value_type> >::result_type
		sqrt  (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_sqrt<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (sin v) [i] = sin (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_sin<typename E::value_type> >::result_type
		sin  (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_sin<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}

	// (cos v) [i] = cos (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits<E, scalar_cos<typename E::value_type> >::result_type
		cos  (const vector_expression<E> &e) {
			typedef typename vector_unary_traits<E, scalar_cos<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (v > t) [i] = v [i] > t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater<typename E1::value_type, typename E1::value_type> >::result_type
		operator > (const vector_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename vector_binary_scalar2_traits<E1, typename E1::value_type, 
				scalar_greater<typename E1::value_type , typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t > v) [i] = t > v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater<typename E2::value_type, typename E2::value_type> >::result_type
		operator > (const typename E2::value_type& e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<typename E2::value_type, E2, 
				scalar_greater<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 > v2) [i] = v1[i] > v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, scalar_greater <typename E1::value_type, typename E2::value_type> >::result_type
		operator > (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, scalar_greater<typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}

	// (v >= t) [i] = v [i] >= t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::result_type
		operator >= (const vector_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename vector_binary_scalar2_traits<E1, typename E1::value_type, 
				scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t >= v) [i] = t >= v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater_equal<typename E2::value_type, typename E2::value_type> >::result_type
		operator >= (const typename E2::value_type& e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<typename E2::value_type, E2, 
				scalar_greater_equal<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 >= v2) [i] = v1[i] >= v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, scalar_greater_equal <typename E1::value_type , typename E2::value_type > >::result_type
		operator >= (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, scalar_greater_equal <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}

	// (v < t) [i] = v [i] < t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser<typename E1::value_type, typename E1::value_type> >::result_type
		operator < (const vector_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename vector_binary_scalar2_traits<E1, typename E1::value_type, 
				scalar_lesser<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t < v) [i] = t < v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser<typename E2::value_type, typename E2::value_type> >::result_type
		operator < (const typename E2::value_type& e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<typename E2::value_type, E2, 
				scalar_lesser<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 < v2) [i] = v1[i] < v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, scalar_lesser <typename E1::value_type, typename E2::value_type> >::result_type
		operator < (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, scalar_lesser<typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (v <= t) [i] = v [i] <= t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::result_type
		operator <= (const vector_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename vector_binary_scalar2_traits<E1, typename E1::value_type, 
				scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t <= v) [i] = t <= v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser_equal<typename E2::value_type, typename E2::value_type> >::result_type
		operator <= (const typename E2::value_type& e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<typename E2::value_type, E2, 
				scalar_lesser_equal<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 <= v2) [i] = v1[i] <= v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, scalar_lesser_equal <typename E1::value_type, typename E2::value_type> >::result_type
		operator <= (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, scalar_lesser_equal<typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (v != t) [i] = v [i] != t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, typename E1::value_type, scalar_inequality<typename E1::value_type, typename E1::value_type> >::result_type
		operator != (const vector_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename vector_binary_scalar2_traits<E1, typename E1::value_type, 
				scalar_inequality<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t != v) [i] = t != v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<typename E2::value_type, E2, scalar_inequality<typename E2::value_type, typename E2::value_type> >::result_type
		operator != (const typename E2::value_type& e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<typename E2::value_type, E2, 
				scalar_inequality<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 != v2) [i] = v1[i] != v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, scalar_inequality <typename E1::value_type, typename E2::value_type> >::result_type
		operator != (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, scalar_inequality<typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (v && t) [i] = v [i] && t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, const bool, logical_and<typename E1::value_type, bool> >::result_type
		operator && (const vector_expression<E1> &e1, const bool &e2) {
			typedef typename vector_binary_scalar2_traits<E1, const bool , logical_and <typename E1::value_type, bool> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t && v) [i] = t && v [i]
	template< class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<const bool, E2, logical_and<bool, typename E2::value_type> >::result_type
		operator && (const bool &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<const bool, E2, logical_and<bool, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 && v2) [i] = v1[i] || v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, logical_and <bool, bool> >::result_type
		operator && (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, logical_and <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}

	// (v || t) [i] = v [i] || t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar2_traits<E1, const bool, logical_or<typename E1::value_type, bool> >::result_type
		operator || (const vector_expression<E1> &e1, const bool &e2) {
			typedef typename vector_binary_scalar2_traits<E1, const bool, logical_or<typename E1::value_type, bool> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t || v) [i] = t || v [i]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_scalar1_traits<const bool, E2, logical_or<bool, typename E2::value_type> >::result_type
		operator || (const bool &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_scalar1_traits<const bool, E2, logical_or<bool, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (v1 || v2) [i] = v1[i] || v2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename vector_binary_traits <E1, E2, logical_or <bool, bool> >::result_type
		operator || (const vector_expression<E1> &e1, const vector_expression<E2> &e2) {
			typedef typename vector_binary_traits<E1, E2, logical_or <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
 	template <class E>
	BOOST_UBLAS_INLINE
		typename vector_unary_traits <E, logical_not < typename E::value_type> > :: result_type
		operator ! (const vector_expression <E> &e) {
			typedef typename vector_unary_traits<E, logical_not <typename E::value_type> > :: expression_type expression_type;
			return expression_type (e () );
	}
}}}
#endif //VECTOR_EXPRESSION_EXT_H
