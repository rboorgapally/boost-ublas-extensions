#ifndef _BOOST_UBLAS_MATRIX_EXPRESSION_EXT_H
#define _BOOST_UBLAS_MATRIX_EXPRESSION_EXT_H

#include "vector_expression_ext.hpp"
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace boost { namespace numeric { namespace ublas {
	template<class E, class F>
	class matrix_vector_unary:
		public vector_expression<matrix_vector_unary<E,F> > {
	public:
		typedef E expression_type;
	private:
		typedef F functor_type;
	public:
		typedef typename E::const_closure_type expression_closure_type;
	private:
		typedef matrix_vector_unary<E,F> self_type;
	public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
		using vector_expression<self_type>::operator ();
#endif
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
		typedef const self_type const_closure_type;
		typedef const_closure_type closure_type;
		typedef unknown_storage_tag storage_category;
		static const unsigned complexity = 1;
		// Construction and destruction
		BOOST_UBLAS_INLINE
			explicit matrix_vector_unary (const expression_type &e):
		e_ (e) {}
		// Accessors
		BOOST_UBLAS_INLINE
			size_type size () const {
				return e_.size2 ();
		}
	public:
		// Expression accessors
		BOOST_UBLAS_INLINE
			const expression_closure_type &expression () const {
				return e_;
		}
	public:
		// Element access
		BOOST_UBLAS_INLINE
			const_reference operator () (size_type i) const {
				return functor_type::apply (e_.size1(), e_.find1(0, 0, i));
		}
		// Closure comparison
		BOOST_UBLAS_INLINE
			bool same_closure (const matrix_vector_unary &mvu) const {
				return (*this).expression ().same_closure (mvu.expression ());
		}
		// Iterator types
	private:
		typedef typename E::const_iterator2 const_subiterator_type;
		typedef const value_type *const_pointer;
	public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
		typedef indexed_const_iterator<const_closure_type, typename const_subiterator_type::iterator_category> const_iterator;
		typedef const_iterator iterator;
#else
		class const_iterator;
		typedef const_iterator iterator;
#endif
		// Element lookup
		BOOST_UBLAS_INLINE
			const_iterator find (size_type i) const {
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
				const_subiterator_type it (e_.find2 (0, 0, i));
				return const_iterator (*this, it.index2 ());
#else
				return const_iterator (*this, e_.find2 (0, 0, i));
#endif
		}
#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
		class const_iterator:
			public container_const_reference<matrix_vector_unary>,
			public iterator_base_traits<typename E::const_iterator2::iterator_category>::template
			iterator_base<const_iterator, value_type>::type {
		public:
			typedef typename E::const_iterator2::iterator_category iterator_category;
			typedef typename matrix_vector_unary::difference_type difference_type;
			typedef typename matrix_vector_unary::value_type value_type;
			typedef typename matrix_vector_unary::const_reference reference;
			typedef typename matrix_vector_unary::const_pointer pointer;
			// Construction and destruction
			BOOST_UBLAS_INLINE
				const_iterator ():
			container_const_reference<self_type> (), it_ () {}
			BOOST_UBLAS_INLINE
				const_iterator (const self_type &mu, const const_subiterator_type &it):
			container_const_reference<self_type> (mu), it_ (it) {}
			// Arithmetic
			BOOST_UBLAS_INLINE
				const_iterator &operator ++ () {
					++ it_;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator &operator -- () {
					-- it_;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator &operator += (difference_type n) {
					it_ += n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator &operator -= (difference_type n) {
					it_ -= n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				difference_type operator - (const const_iterator &it) const {
					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
					return it_ - it.it_;
			}
			// Dereference
			BOOST_UBLAS_INLINE
				const_reference operator * () const {
					return functor_type::apply ((*this) ().expression ().size1 (), it_.begin());
			}
			// Indices
			BOOST_UBLAS_INLINE
				size_type index () const {
					return it_.index();
			}
			// Assignment 
			BOOST_UBLAS_INLINE
				const_iterator &operator = (const const_iterator &it) {
					container_const_reference<self_type>::assign (&it ());
					it_ = it.it_;
					return *this;
			}
			// Comparison
			BOOST_UBLAS_INLINE
				bool operator == (const const_iterator &it) const {
					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
					return it_ == it.it_;
			}
			BOOST_UBLAS_INLINE
				bool operator < (const const_iterator &it) const {
					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
					return it_ < it.it_;
			}
		private:
			const_subiterator_type it_;
		};
#endif
		BOOST_UBLAS_INLINE
			const_iterator begin () const {
				return find (0);
		}
		BOOST_UBLAS_INLINE
			const_iterator end () const {
				return find (size ()); 
		}
		// Reverse iterator
		typedef reverse_iterator_base<const_iterator> const_reverse_iterator;
		BOOST_UBLAS_INLINE
			const_reverse_iterator rbegin () const {
				return const_reverse_iterator (end ());
		}
		BOOST_UBLAS_INLINE
			const_reverse_iterator rend () const {
				return const_reverse_iterator (begin ());
		}
	private:
		expression_closure_type e_;
	};
	template <class E1, class F>
	struct matrix_vector_unary_traits {
		typedef unknown_storage_tag storage_category;
		typedef row_major_tag orientation_category;
		typedef typename E1::value_type value_type;
		//typedef typename promote_traits<value_type, value_type>::promote_type promote_type;
		typedef matrix_vector_unary<E1, F> expression_type;
#ifndef BOOST_UBLAS_SIMPLE_ET_DEBUG
		typedef expression_type result_type;
#else
		typedef typename E1::vector_temporary_type result_type;
#endif
	};
	// Matrix returning vector
	template<class E>
	struct matrix_vector_unary_functor {
		typedef std::size_t size_type;
		typedef std::ptrdiff_t difference_type;
		typedef typename E::value_type value_type;
		typedef typename E::value_type result_type;
		typedef typename E::const_iterator1 const_iterator1;
	};
	template<class E>
	struct matrix_column_sum: 
		public matrix_vector_unary_functor<E> {
			typedef typename matrix_vector_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size, const const_iterator1 &it) { 
					result_type t = result_type (0);
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i)
						t += *it_temp, ++ it_temp;
					return t;
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) { 
					result_type t = result_type (0);
					while (-- size >= 0)
						t += *it, ++ it;
					return t; 
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) {
					result_type t = result_type (0);
					while (it != it_end) 
						t += *it, ++ it;
					return t; 
			}
	};
	//sum m[j] = sum (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_column_sum< E> >::result_type
		sum (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_column_sum<E> >::expression_type expression_type;
			return expression_type (e ());
	}
	//norm_L1 m[j] = norm_L1 (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_norm_1 <E> >::result_type
		norm_L1 (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_1 <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	//norm_L2 m[j] = norm_L2 (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_norm_2 <E> >::result_type
		norm_L2 (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_2 <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	//norm_inf m[j] = norm_inf (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_norm_inf <E> >::result_type
		col_norm_inf (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_inf <E> >::expression_type expression_type;
			return expression_type (e ());
	}

	//max m[j] = max (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_max <E> >::result_type
		max (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_max <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	struct matrix_vector_min:
		public matrix_vector_real_unary_functor<E> {
			typedef typename matrix_vector_real_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<E>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size,const const_iterator1 &it) {
					real_type t = (std::numeric_limits<result_type>::max());
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i) {
						real_type u = *it_temp;
						if (u < t)
							t = u;
						++it_temp;		
					}
					return t;
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = (std::numeric_limits<result_type>::max());
					while (-- size >= 0) {
						real_type u = *it;
						if (u < t)
							t = u;
						++ it;
					}
					return t;
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = (std::numeric_limits<result_type>::max());
					while (it != it_end) {
						real_type u = *it;
						if (u < t) 
							t = u;
						++ it;
					}
					return t; 
			}
	};
	//min m[j] = min (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_min <E> >::result_type
		min (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_min <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	struct matrix_vector_abs_min:
		public matrix_vector_real_unary_functor<E> {
			typedef typename matrix_vector_real_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<E>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size,const const_iterator1 &it) {
					real_type t = (std::numeric_limits<result_type>::max());
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i) {
						real_type u = std::abs(*it_temp);
						if (u < t)
							t = u;
						++it_temp;		
					}
					return t;
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = (std::numeric_limits<result_type>::max());
					while (-- size >= 0) {
						real_type u = std::abs(*it);
						if (u < t)
							t = u;
						++ it;
					}
					return t;
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = (std::numeric_limits<result_type>::max());
					while (it != it_end) {
						real_type u = std::abs(*it);
						if (u < t) 
							t = u;
						++ it;
					}
					return t; 
			}
	};
	//abs_min m[j] = abs_min (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_abs_min <E> >::result_type
		abs_min (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_abs_min <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	struct matrix_vector_abs_max:
		public matrix_vector_real_unary_functor<E> {
			typedef typename matrix_vector_real_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<E>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size,const const_iterator1 &it) {
					real_type t = -(std::numeric_limits<result_type>::max());
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i) {
						real_type u = std::abs(*it_temp);
						if (u > t)
							t = u;
						++it_temp;		
					}
					return t;
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = real_type ();
					while (-- size >= 0) {
						real_type u = std::abs(*it);
						if (u > t)
							t = u;
						++ it;
					}
					return t;
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = real_type ();
					while (it != it_end) {
						real_type u = std::abs(*it);
						if (u > t) 
							t = u;
						++ it;
					}
					return t; 
			}
	};
	//abs_max m[j] = abs_max (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_abs_max <E> >::result_type
		abs_max (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_abs_max <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	struct matrix_vector_mean:
		public matrix_vector_real_unary_functor<E> {
			typedef typename matrix_vector_real_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<E>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size, const const_iterator1 &it) {
					result_type t = result_type (0);
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i)
						t += *it_temp, ++ it_temp;
					return (t/size);
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = real_type (0);
					while (-- size >= 0) {
						real_type u = *it;
						t += u;
						++ it;
					}
					return (t/size);
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = real_type (0);
					int size = 0;
					while (it != it_end) {
						real_type u = *it;
						t += u;
						++ it;
						++ size;
					}
					return (t/size); 
			}
	};
	//mean m[j] = mean (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_mean <E> >::result_type
		mean (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_mean <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	template<class E>
	struct matrix_vector_abs_mean:
		public matrix_vector_real_unary_functor<E> {
			typedef typename matrix_vector_real_unary_functor<E>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<E>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<E>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<E>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<E>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<E>::const_iterator1 const_iterator1;
			//template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size, const const_iterator1 &it) {
					result_type t = result_type (0);
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i)
						t += std::abs(*it_temp), ++ it_temp;
					return (t/size);
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = real_type ();
					while (-- size >= 0) {
						real_type u = std::abs(*it);
						t += u;
						++ it;
					}
					return (t/size);
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = real_type ();
					int size = 0;
					while (it != it_end) {
						real_type u = std::abs(*it);
						t += u;
						++ it;
						++ size;
					}
					return (t/size); 
			}
	};
	//mean m[j] = mean (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_abs_mean <E> >::result_type
		abs_mean (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_abs_mean <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	//rms m[j] = rms (m[i][j])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_vector_unary_traits<E, matrix_vector_rms <E> >::result_type
		rms (const matrix_expression<E> &e) {
			typedef typename matrix_vector_unary_traits<E, matrix_vector_rms <E> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (t + m) [i] [j] = t + m [i] [j]
	template<class T1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::result_type
		operator + (const T1 &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<const T1, E2, scalar_plus<T1, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m + t) [i] [j] = m [i] [j] + t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::result_type
		operator + (const matrix_expression<E1> &e1, const T2 &e2) {
			typedef typename matrix_binary_scalar2_traits<E1, const T2, scalar_plus<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1(), e2);
	}
	// (t - m) [i] [j] = t - m [i] [j]
	template<class T1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<const T1, E2, scalar_minus<T1, typename E2::value_type> >::result_type
		operator - (const T1 &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<const T1, E2, scalar_minus<T1, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m - t) [i] [j] = m [i] [j] - t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, const T2, scalar_minus<typename E1::value_type, T2> >::result_type
		operator - (const matrix_expression<E1> &e1, const T2 &e2) {
			typedef typename matrix_binary_scalar2_traits<E1, const T2, scalar_minus<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1(), e2);
	}

	// (m ^ t) [i] [j] = m [i] [j] ^ t
	template<class E1, class T2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::result_type
		operator ^ (const matrix_expression<E1> &e1, const T2 &e2) {
			typedef typename matrix_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::expression_type expression_type;
			return expression_type (e1(), e2);
	}
	// (log10 v) [i] = log10 (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_log10<typename E::value_type> >::result_type
		log10 (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_log10<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (abs v) [i] = abs (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_abs<typename E::value_type> >::result_type
		abs (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_abs<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (exp v) [i] = exp (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::result_type
		exp (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (log v) [i] = log (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::result_type
		log (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (sqrt v) [i] = sqrt (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_sqrt<typename E::value_type> >::result_type
		sqrt  (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_sqrt<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (sin v) [i] = sin (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_sin<typename E::value_type> >::result_type
		sin  (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_sin<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (cos v) [i] = cos (v [i])
	template<class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits<E, scalar_cos<typename E::value_type> >::result_type
		cos  (const matrix_expression<E> &e) {
			typedef typename matrix_unary1_traits<E, scalar_cos<typename E::value_type> >::expression_type expression_type;
			return expression_type (e ());
	}
	// (t && m) [i] [j] = t && m [i] [j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits< bool, E2, logical_and<bool, typename E2::value_type> >::result_type
		operator && (bool e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits< bool, E2, logical_and<bool, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m && t) [i] [j] = m [i] [j] && t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, bool, logical_and<typename E1::value_type, bool> >::result_type
		operator && (const matrix_expression<E1> &e1,  bool e2) {
			typedef typename matrix_binary_scalar2_traits<E1, bool, logical_and<typename E1::value_type, bool> >::expression_type expression_type;
			return expression_type (e1(), e2);
	}
	//(m1 && m2) [i] [j] = m1 [i] [j] && m2 [i] [j]
	template<class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits<E1, E2, logical_and<typename E1::value_type, typename E2::value_type> >::result_type
		operator && (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, logical_and<typename E1::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1(), e2());
	}
	// (t || m) [i] [j] = t || m [i] [j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<bool, E2, logical_or<bool, typename E2::value_type> >::result_type
		operator || (bool e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<bool, E2, logical_or<bool, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m || t) [i] [j] = m [i] [j] || t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, bool, logical_or<typename E1::value_type, bool> >::result_type
		operator || (const matrix_expression<E1> &e1, bool e2) {
			typedef typename matrix_binary_scalar2_traits<E1, bool, logical_or<typename E1::value_type, bool> >::expression_type expression_type;
			return expression_type (e1(), e2);
	}
	//(m1 || m2) [i] [j] = m1 [i] [j] || m2 [i] [j]
	template<class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits<E1, E2, logical_or<typename E1::value_type, typename E2::value_type> >::result_type
		operator || (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, logical_or<typename E1::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1(), e2());
	}
	template <class E>
	BOOST_UBLAS_INLINE
		typename matrix_unary1_traits <E, logical_not < typename E::value_type> > :: result_type
		operator ! (const matrix_expression <E> &e) {
			typedef typename matrix_unary1_traits<E, logical_not <typename E::value_type> > :: expression_type expression_type;
			return expression_type (e () );
	}
	// ((m > t) [i][j]  =  m [i][j] > t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1,typename E1::value_type, scalar_greater<typename E1::value_type, typename E1::value_type> >::result_type
		operator > (const matrix_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t > m) [i][j]  = t > m [i][j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater<typename E2::value_type, typename E2::value_type> >::result_type
		operator > (const typename E2::value_type& e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m1 > m2) [i][j] = m1[i][j] > m2[i][j]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits <E1, E2, scalar_greater <typename E1::value_type, typename E2::value_type> >::result_type
		operator > (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, scalar_greater<typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (m >= t) [i][j] = m [i][j] >= t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::result_type
		operator >= (const matrix_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t >= m) [i][j] = t >= v [i][j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater_equal<typename E2::value_type, typename E2::value_type> >::result_type
		operator >= (const typename E2::value_type& e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater_equal<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m1 >= m2) [i] = m1[i] >= m2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits <E1, E2, scalar_greater_equal <typename E1::value_type , typename E2::value_type > >::result_type
		operator >= (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, scalar_greater_equal <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (m < t) [i][j] = m [i][j] < t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser<typename E1::value_type, typename E1::value_type> >::result_type
		operator < (const matrix_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (m < v) [i][j] = t < m [i][j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser<typename E2::value_type, typename E2::value_type> >::result_type
		operator < (const typename E2::value_type& e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m1 < m2) [i] = m1[i] < m2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits <E1, E2, scalar_lesser <typename E1::value_type , typename E2::value_type > >::result_type
		operator < (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, scalar_lesser <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (m <= t) [i][j] = m [i][j] <= t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::result_type
		operator <= (const matrix_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t <= m) [i][j] = t <= m [i][j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser_equal<typename E2::value_type, typename E2::value_type> >::result_type
		operator <= (const typename E2::value_type& e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser_equal<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m1 <= m2) [i] = m1[i] <= m2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits <E1, E2, scalar_lesser_equal <typename E1::value_type , typename E2::value_type > >::result_type
		operator <= (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, scalar_lesser_equal <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
	// (m != t) [i][j] = m [i][j] != t
	template<class E1>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_inequality<typename E1::value_type, typename E1::value_type> >::result_type
		operator != (const matrix_expression<E1> &e1, const typename E1::value_type& e2) {
			typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_inequality<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
			return expression_type (e1 (), e2);
	}
	// (t != m) [i][j] = t != m [i][j]
	template<class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_inequality<typename E2::value_type, typename E2::value_type> >::result_type
		operator != (const typename E2::value_type& e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_inequality<typename E2::value_type, typename E2::value_type> >::expression_type expression_type;
			return expression_type (e1, e2 ());
	}
	// (m1 != m2) [i] = m1[i] != m2[i]
	template <class E1, class E2>
	BOOST_UBLAS_INLINE
		typename matrix_binary_traits <E1, E2, scalar_inequality <typename E1::value_type , typename E2::value_type > >::result_type
		operator != (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2) {
			typedef typename matrix_binary_traits<E1, E2, scalar_inequality <typename E1::value_type, typename E2::value_type> >:: expression_type expression_type;
			return expression_type (e1 (), e2 ());
	}
//	template<class E1, class E2>
//	class vector_matrix_binary_cat:
//		public matrix_expression<vector_matrix_binary_cat<E1, E2, F> > {
//			typedef E1 expression1_type;
//			typedef E2 expression2_type;
//	public:
//		typedef typename E1::const_closure_type expression1_closure_type;
//		typedef typename E2::const_closure_type expression2_closure_type;
//	private:
//		typedef vector_matrix_binary_cat<E1, E2> self_type;
//	public:
//#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
//		using matrix_expression<self_type>::operator ();
//#endif
//		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
//		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
//		typedef typename promote_traits<typename E1::value_type, typename E2::value_type>::promote_type value_type;
//		typedef value_type const_reference;
//		typedef const_reference reference;
//		typedef const self_type const_closure_type;
//		typedef const_closure_type closure_type;
//		typedef unknown_orientation_tag orientation_category;
//		typedef unknown_storage_tag storage_category;
//		// Construction and destruction 
//		BOOST_UBLAS_INLINE
//			vector_matrix_binary_cat (const expression1_type &e1, const expression2_type &e2): 
//		e1_ (e1), e2_ (e2) {}
//		// Accessors
//		BOOST_UBLAS_INLINE
//			size_type size1 () const {
//				BOOST_UBLAS_SAME(e1_.size(), e2_.size1());
//		}
//		BOOST_UBLAS_INLINE
//			size_type size2 () const { 
//				return e2_.size2 ();
//		}
//	public:
//		// Expression accessors
//		BOOST_UBLAS_INLINE
//			const expression1_closure_type &expression1 () const {
//				return e1_;
//		}
//		BOOST_UBLAS_INLINE
//			const expression2_closure_type &expression2 () const {
//				return e2_;
//		}
//	public:
//		// Element access
//		BOOST_UBLAS_INLINE
//			const_reference operator () (size_type i, size_type j) const {
//				return (j > 0)?e2_(i, j -1):e1_(i);
//		}
//		// Closure comparison
//		BOOST_UBLAS_INLINE
//			bool same_closure (const vector_matrix_binary_cat &vmb) const {
//				return (*this).expression1 ().same_closure (vmb.expression1 ()) &&
//					(*this).expression2 ().same_closure (vmb.expression2 ());
//		}
//		// Iterator types
//	private:
//		typedef typename E1::const_iterator const_subiterator1_type;
//		typedef typename E2::const_iterator const_subiterator2_type;
//		typedef const value_type *const_pointer;
//	public:
//#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
//		typedef typename iterator_restrict_traits<typename const_subiterator1_type::iterator_category,
//			typename const_subiterator2_type::iterator_category>::iterator_category iterator_category;
//		typedef indexed_const_iterator1<const_closure_type, iterator_category> const_iterator1;
//		typedef const_iterator1 iterator1;
//		typedef indexed_const_iterator2<const_closure_type, iterator_category> const_iterator2;
//		typedef const_iterator2 iterator2;
//#else
//		class const_iterator1;
//		typedef const_iterator1 iterator1;
//		class const_iterator2;
//		typedef const_iterator2 iterator2;
//#endif
//		typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
//		typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
//		// Element lookup
//		BOOST_UBLAS_INLINE
//			const_iterator1 find1 (int rank, size_type i, size_type j) const {
//				const_subiterator1_type it1 (e1_.find (i));
//				const_subiterator1_type it1_end (e1_.find (size1 ()));
//				const_subiterator2_type it2 (e2_.find (j));
//				const_subiterator2_type it2_end (e2_.find (size2 ()));
//				if (it2 == it2_end || (rank == 1 && (it2.index () != j || *it2 == value_type/*zero*/()))) {
//					it1 = it1_end;
//					it2 = it2_end;
//				}
//#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
//				return const_iterator1 (*this, it1.index (), it2.index ());
//#else
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//				return const_iterator1 (*this, it1, it2, it2 != it2_end ? *it2 : value_type/*zero*/());
//#else
//				return const_iterator1 (*this, it1, it2);
//#endif
//#endif
//		}
//		BOOST_UBLAS_INLINE
//			const_iterator2 find2 (int rank, size_type i, size_type j) const {
//				const_subiterator2_type it2 (e2_.find (j));
//				const_subiterator2_type it2_end (e2_.find (size2 ()));
//				const_subiterator1_type it1 (e1_.find (i));
//				const_subiterator1_type it1_end (e1_.find (size1 ()));
//				if (it1 == it1_end || (rank == 1 && (it1.index () != i || *it1 == value_type/*zero*/()))) {
//					it2 = it2_end;
//					it1 = it1_end;
//				}
//#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
//				return const_iterator2 (*this, it1.index (), it2.index ());
//#else
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//				return const_iterator2 (*this, it1, it2, it1 != it1_end ? *it1 : value_type/*zero*/());
//#else
//				return const_iterator2 (*this, it1, it2);
//#endif
//#endif
//		}
//		// Iterators enhance the iterators of the referenced expressions
//		// with the binary functor.
//#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
//		class const_iterator1:
//			public container_const_reference<vector_matrix_binary_cat>,
//			public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
//			typename E2::const_iterator::iterator_category>::iterator_category>::template
//			iterator_base<const_iterator1, value_type>::type {
//		public:
//			typedef typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
//				typename E2::const_iterator::iterator_category>::iterator_category iterator_category;
//			typedef typename vector_matrix_binary_cat::difference_type difference_type;
//			typedef typename vector_matrix_binary_cat::value_type value_type;
//			typedef typename vector_matrix_binary_cat::const_reference reference;
//			typedef typename vector_matrix_binary_cat::const_pointer pointer;
//			typedef const_iterator2 dual_iterator_type;
//			typedef const_reverse_iterator2 dual_reverse_iterator_type;
//			// Construction and destruction
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//			BOOST_UBLAS_INLINE
//				const_iterator1 ():
//			container_const_reference<self_type> (), it1_ (), it2_ (), t2_ () {}
//			BOOST_UBLAS_INLINE
//				const_iterator1 (const self_type &vmb, const const_subiterator1_type &it1, const const_subiterator2_type &it2, value_type t2):
//			container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2), t2_ (t2) {}
//#else
//			BOOST_UBLAS_INLINE
//				const_iterator1 ():
//			container_const_reference<self_type> (), it1_ (), it2_ () {}
//			BOOST_UBLAS_INLINE
//				const_iterator1 (const self_type &vmb, const const_subiterator1_type &it1, const const_subiterator2_type &it2):
//			container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2) {}
//#endif
//			// Arithmetic
//			BOOST_UBLAS_INLINE
//				const_iterator1 &operator ++ () {
//					++ it1_;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator1 &operator -- () {
//					-- it1_;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator1 &operator += (difference_type n) {
//					it1_ += n;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator1 &operator -= (difference_type n) {
//					it1_ -= n;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				difference_type operator - (const const_iterator1 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
//					return it1_ - it.it1_;
//			}
//			// Dereference
//			BOOST_UBLAS_INLINE
//				const_reference operator * () const {
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//					return functor_type::apply (*it1_, t2_);
//#else
//					return functor_type::apply (*it1_, *it2_);
//#endif
//			}
//			BOOST_UBLAS_INLINE
//				const_reference operator [] (difference_type n) const {
//					return *(*this + n);
//			}
//#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_iterator2 begin () const {
//					return (*this) ().find2 (1, index1 (), 0);
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_iterator2 end () const {
//					return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_reverse_iterator2 rbegin () const {
//					return const_reverse_iterator2 (end ());
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_reverse_iterator2 rend () const {
//					return const_reverse_iterator2 (begin ());
//			}
//#endif
//			// Indices
//			BOOST_UBLAS_INLINE
//				size_type index1 () const {
//					return it1_.index ();
//			}
//			BOOST_UBLAS_INLINE
//				size_type  index2 () const {
//					return it2_.index ();
//			}
//			// Assignment
//			BOOST_UBLAS_INLINE
//				const_iterator1 &operator = (const const_iterator1 &it) {
//					container_const_reference<self_type>::assign (&it ());
//					it1_ = it.it1_;
//					it2_ = it.it2_;
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//					t2_ = it.t2_;
//#endif
//					return *this;
//			}
//			// Comparison
//			BOOST_UBLAS_INLINE
//				bool operator == (const const_iterator1 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
//					return it1_ == it.it1_;
//			}
//			BOOST_UBLAS_INLINE
//				bool operator < (const const_iterator1 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it2_ == it.it2_, external_logic ());
//					return it1_ < it.it1_;
//			}
//		private:
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//			const_subiterator1_type it1_;
//			// Mutable due to assignment
//			/* const */ const_subiterator2_type it2_;
//			value_type t2_;
//#else
//			const_subiterator1_type it1_;
//			const_subiterator2_type it2_;
//#endif
//		};
//#endif
//		BOOST_UBLAS_INLINE
//			const_iterator1 begin1 () const {
//				return find1 (0, 0, 0);
//		}
//		BOOST_UBLAS_INLINE
//			const_iterator1 end1 () const {
//				return find1 (0, size1 (), 0);
//		}
//#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
//		class const_iterator2:
//			public container_const_reference<vector_matrix_binary_cat>,
//			public iterator_base_traits<typename iterator_restrict_traits<typename E1::const_iterator::iterator_category,
//			typename E2::const_iterator::iterator_category>::iterator_category>::template
//			iterator_base<const_iterator2, value_type>::type {
//		public:
//			typedef typename iterator_restrict_traits<typename E1::const_iterator::iterator_category, 
//				typename E2::const_iterator::iterator_category>::iterator_category iterator_category;
//			typedef typename vector_matrix_binary_cat::difference_type difference_type;
//			typedef typename vector_matrix_binary_cat::value_type value_type;
//			typedef typename vector_matrix_binary_cat::const_reference reference;
//			typedef typename vector_matrix_binary_cat::const_pointer pointer;
//			typedef const_iterator1 dual_iterator_type;
//			typedef const_reverse_iterator1 dual_reverse_iterator_type;
//			// Construction and destruction
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//			BOOST_UBLAS_INLINE
//				const_iterator2 ():
//			container_const_reference<self_type> (), it1_ (), it2_ (), t1_ () {}
//			BOOST_UBLAS_INLINE
//				const_iterator2 (const self_type &vmb, const const_subiterator1_type &it1, const const_subiterator2_type &it2, value_type t1):
//			container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2), t1_ (t1) {}
//#else
//			BOOST_UBLAS_INLINE
//				const_iterator2 ():
//			container_const_reference<self_type> (), it1_ (), it2_ () {}
//			BOOST_UBLAS_INLINE
//				const_iterator2 (const self_type &vmb, const const_subiterator1_type &it1, const const_subiterator2_type &it2):
//			container_const_reference<self_type> (vmb), it1_ (it1), it2_ (it2) {}
//#endif
//			// Arithmetic
//			BOOST_UBLAS_INLINE
//				const_iterator2 &operator ++ () {
//					++ it2_;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator2 &operator -- () {
//					-- it2_;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator2 &operator += (difference_type n) {
//					it2_ += n;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				const_iterator2 &operator -= (difference_type n) {
//					it2_ -= n;
//					return *this;
//			}
//			BOOST_UBLAS_INLINE
//				difference_type operator - (const const_iterator2 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure(it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
//					return it2_ - it.it2_;
//			}
//			// Dereference
//			BOOST_UBLAS_INLINE
//				const_reference operator * () const {
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//					return functor_type::apply (t1_, *it2_);
//#else
//					return functor_type::apply (*it1_, *it2_);
//#endif
//			}
//			BOOST_UBLAS_INLINE
//				const_reference operator [] (difference_type n) const {
//					return *(*this + n);
//			}
//#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_iterator1 begin () const {
//					return (*this) ().find1 (1, 0, index2 ());
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_iterator1 end () const {
//					return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_reverse_iterator1 rbegin () const {
//					return const_reverse_iterator1 (end ());
//			}
//			BOOST_UBLAS_INLINE
//#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
//				typename self_type::
//#endif
//				const_reverse_iterator1 rend () const {
//					return const_reverse_iterator1 (begin ());
//			}
//#endif
//			// Indices
//			BOOST_UBLAS_INLINE
//				size_type index1 () const {
//					return it1_.index ();
//			}
//			BOOST_UBLAS_INLINE
//				size_type  index2 () const {
//					return it2_.index ();
//			}
//			// Assignment
//			BOOST_UBLAS_INLINE
//				const_iterator2 &operator = (const const_iterator2 &it) {
//					container_const_reference<self_type>::assign (&it ());
//					it1_ = it.it1_;
//					it2_ = it.it2_;
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//					t1_ = it.t1_;
//#endif
//					return *this;
//			}
//			// Comparison
//			BOOST_UBLAS_INLINE
//				bool operator == (const const_iterator2 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure( it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
//					return it2_ == it.it2_;
//			}
//			BOOST_UBLAS_INLINE
//				bool operator < (const const_iterator2 &it) const {
//					BOOST_UBLAS_CHECK ((*this) ().same_closure (it ()), external_logic ());
//					BOOST_UBLAS_CHECK (it1_ == it.it1_, external_logic ());
//					return it2_ < it.it2_;
//			}
//		private:
//#ifdef BOOST_UBLAS_USE_INVARIANT_HOISTING
//			// Mutable due to assignment
//			/* const */ const_subiterator1_type it1_;
//			const_subiterator2_type it2_;
//			value_type t1_;
//#else
//			const_subiterator1_type it1_;
//			const_subiterator2_type it2_;
//#endif
//		};
//#endif
//		BOOST_UBLAS_INLINE
//			const_iterator2 begin2 () const {
//				return find2 (0, 0, 0);
//		}
//		BOOST_UBLAS_INLINE
//			const_iterator2 end2 () const {
//				return find2 (0, 0, size2 ());
//		}
//		// Reverse iterators
//		BOOST_UBLAS_INLINE
//			const_reverse_iterator1 rbegin1 () const {
//				return const_reverse_iterator1 (end1 ());
//		}
//		BOOST_UBLAS_INLINE
//			const_reverse_iterator1 rend1 () const {
//				return const_reverse_iterator1 (begin1 ());
//		}
//		BOOST_UBLAS_INLINE
//			const_reverse_iterator2 rbegin2 () const {
//				return const_reverse_iterator2 (end2 ());
//		}
//		BOOST_UBLAS_INLINE
//			const_reverse_iterator2 rend2 () const {
//				return const_reverse_iterator2 (begin2 ());
//		}
//	private:
//		expression1_closure_type e1_;
//		expression2_closure_type e2_;
//	};
//	template<class E1, class E2, class F>
//	struct vector_matrix_binary_cat_traits {
//		typedef vector_matrix_binary_cat<E1, E2, F> expression_type;
//#ifndef BOOST_UBLAS_SIMPLE_ET_DEBUG
//		typedef expression_type result_type; 
//#else
//		// ISSUE matrix is arbitary temporary type
//		typedef matrix<typename F::value_type> result_type;
//#endif
//	};
}}}
#endif  //_BOOST_UBLAS_MATRIX_EXPRESSION_EXT_H
