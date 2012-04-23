#ifndef FUNCTIONAL_EXT_H
#define FUNCTIONAL_EXT_H

#include <boost/numeric/ublas/functional.hpp>
namespace boost { namespace numeric { namespace ublas {

	template<class T, class OP>
	struct general_vector_scalar_unary: 
		public vector_scalar_unary_functor<T> {

			//typedef typename vector_scalar_unary_functor<T>::difference_type difference_type;
			typedef typename vector_scalar_unary_functor<T>::value_type value_type;
			typedef typename vector_scalar_unary_functor<T>::result_type result_type;

			// general case, access by index
			template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (const vector_expression<E> &e) { 
					result_type t = OP::initial_value();
					typedef typename E::size_type size_type;
					size_type size (e ().size ());
					for (size_type i = 0; i < size; ++ i)
						OP::update(t, e () (i));
					return t;
			};
			// Dense case
			template<class D,class I>
			static BOOST_UBLAS_INLINE
				result_type apply (D size, I it) { 
					result_type t = OP::initial_value();
					while (-- size >= 0)
						OP::update(t, *it), ++ it;
					return t; 
			}

			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) {
					result_type t = OP::initial_value();
					while (it != it_end) 
						OP::update(t, *it), ++ it;
					return t; 
			}
	};


	template <class T>
	struct scalar_log_minimum {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return std::log(std::numeric_limits<result_type>::max());
		}
		static
			void update(result_type& t, const value_type& x){
				if ( (x>0) && (std::log(x)<t) ) t = std::log(x);
		}
	};
	template <class T>
	struct scalar_log_maximum {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return T(1.0);
		}
		static
			void update(result_type& t, const value_type& x){
				if ( (x>0) && (std::log(x)>t) ) t = std::log(x);
		}
	};
	template <class T> 
	struct scalar_minimum { 
	public: 
		typedef T value_type; 
		typedef T result_type; 
		static 
			result_type initial_value(){
				//std::cout << std::numeric_limits<result_type>::max() << std::endl;
				return (std::numeric_limits<result_type>::max()); 
		}
		static 
			void update(result_type& t, const value_type& x){ 
				if (x<t) t=x; 
		} 
	}; 
	template <class T>
	struct scalar_maximum {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				//std::cout << std::numeric_limits<result_type>::min() << std::endl;
				return -(std::numeric_limits<result_type>::max());
		}
		static
			void update(result_type& t, const value_type& x){
				if (x>t) t=x;
		}
	};
	//The minimum value should be zero
	template <class T>
	struct sum_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static 
			result_type initial_value(){
				return 0;
		}
		static
			void update(result_type& t, const value_type& x){
				t= t+x; 
		}
	};
	template <class T>
	struct abs_min_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return (std::numeric_limits<result_type>::max());
		}
		static
			void update(result_type& t, const value_type& x){
				if (std::abs(x) < t) t = std::abs(x);
		}
	};
	template <class T>
	struct abs_max_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return -(std::numeric_limits<result_type>::max());
		}
		static
			void update(result_type& t, const value_type& x){
				if (std::abs(x) > t) t = std::abs(x);
		}
	};
	template <class T>
	struct mean_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return 0;
		}
		static
			void update(result_type& t, const value_type& x){
				static int num = 0;
				t = (((num*t) + x)/(num+1));
				num++;
		}
	};
	template <class T>
	struct abs_mean_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return 0;
		}
		static
			void update(result_type& t, const value_type& x){
				static int num = 0;
				t = (((num*t) + std::abs(x))/(num+1));
				num++;
		}
	};
	template <class T>
	struct rms_vector {
	public:
		typedef T value_type;
		typedef T result_type;
		static
			result_type initial_value(){
				return 0;
		}
		static
			void update(result_type& t, const value_type& x){
				t = std::sqrt(pow(t,(int)2)+ pow(x,(int)2)); 
		}
	};


	template<class T1, class T2>
	struct scalar_power:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return std::pow(t1,t2);
			}
	};

	template<class T>
	struct scalar_log10:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::log10(t);
			}
	};
	template<class T>
	struct scalar_abs:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::abs(t);
			}
	};
	template<class T>
	struct scalar_exp:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::exp(t);
			}
	};
	template<class T>
	struct scalar_log:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					if (t > 0)return std::log(t); // insert any function you want
					else return 0;			
			}

	};
	template<class T>
	struct scalar_sqrt:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::sqrt(t);
			}
	};
	template<class T>
	struct scalar_sin:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::sin (t);
			}
	};
	template<class T>
	struct scalar_cos:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			typedef typename scalar_unary_functor<T>::result_type result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return std::cos (t);
			}
	};
	template<class T1, class T2>
	struct logical_and:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 && t2);
			}
	};
	template<class T1, class T2>
	struct logical_or:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 || t2);
			}
	};
	template<class T>
	struct logical_not:
		public scalar_unary_functor<T> {
			typedef typename scalar_unary_functor<T>::argument_type argument_type;
			//typedef typename scalar_unary_functor<T>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument_type t) {
					return (!t);
			}
	};
	template<class T1, class T2>
	struct scalar_greater:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 > t2);
			}
	};
	template<class T1, class T2>
	struct scalar_greater_equal:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 >= t2);
			}
	};
	template<class T1, class T2>
	struct scalar_lesser:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 < t2);
			}
	};
	template<class T1, class T2>
	struct scalar_lesser_equal:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;
			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 <= t2);
			}
	};

	template<class T1, class T2>
	struct scalar_inequality:
		public scalar_binary_functor<T1, T2> {
			typedef typename scalar_binary_functor<T1, T2>::argument1_type argument1_type;
			typedef typename scalar_binary_functor<T1, T2>::argument2_type argument2_type;
			//typedef typename scalar_binary_functor<T1, T2>::result_type result_type;
			typedef bool result_type;

			static BOOST_UBLAS_INLINE
				result_type apply (argument1_type t1, argument2_type t2) {
					return (t1 != t2);
			}
	};
	// Unary returning real scalar 
	template<class E>
	struct matrix_vector_real_unary_functor {
		typedef std::size_t size_type;
		typedef std::ptrdiff_t difference_type;
		typedef typename E::value_type value_type;
		typedef typename type_traits<typename E::value_type>::real_type real_type;
		typedef real_type result_type;
		typedef typename E::const_iterator1 const_iterator1;
	};
	template<class E>
	struct matrix_vector_norm_1:
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
					real_type t = real_type ();
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i) {
						real_type u (std::abs (*it_temp));
						t += u;
						++ it_temp;
					}
					return t;
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = real_type ();
					while (-- size >= 0) {
						real_type u (std::abs(*it));
						t += u;
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
						real_type u (std::abs (*it));
						t += u;
						++ it;
					}
					return t;
			}
	};
	template<class E>
	struct matrix_vector_norm_2:
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
#ifndef BOOST_UBLAS_SCALED_NORM
					real_type t = real_type ();
					const_iterator1 it_temp = it;
					for (size_type i = 0; i < size; ++ i) {
						real_type u (std::abs (*it_temp));
						t +=  u * u, ++ it_temp;
					}
					return std::sqrt (t);
#else
					real_type scale = real_type ();
					real_type sum_squares (1);
					for (size_type i = 0; i < size; ++ i) {
						real_type u (std::abs (*it));
						if (scale < u) {
							real_type v (scale / u);
							sum_squares = sum_squares * v * v + real_type (1);
							scale = u;
						} else {
							real_type v (u / scale);
							sum_squares += v * v;
						}
					}
					return scale * std::sqrt (sum_squares);
#endif
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
#ifndef BOOST_UBLAS_SCALED_NORM
					real_type t = real_type ();
					while (-- size >= 0) {
						real_type u (std::abs (*it));
						t +=  u * u;
						++ it;
					}
					return std::sqrt (t);
#else
					real_type scale = real_type ();
					real_type sum_squares (1);
					while (-- size >= 0) {
						real_type u (std::abs (*it));
						if (scale < u) {
							real_type v (scale / u);
							sum_squares = sum_squares * v * v + real_type (1);
							scale = u;
						} else {
							real_type v (u / scale);
							sum_squares += v * v;
						}
						++ it;
					}
					return scale * std::sqrt (sum_squares);
#endif
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) {
#ifndef BOOST_UBLAS_SCALED_NORM
					real_type t = real_type ();
					while (it != it_end) {
						real_type u (std::abs (*it));
						t +=  u * u;
						++ it;
					}
					return std::sqrt (t);
#else
					real_type scale = real_type ();
					real_type sum_squares (1);
					while (it != it_end) {
						real_type u (std::abs (*it));
						if (scale < u) {
							real_type v (scale / u);
							sum_squares = sum_squares * v * v + real_type (1);
							scale = u;
						} else {
							real_type v (u / scale);
							sum_squares += v * v;
						}
						++ it;
					}
					return scale * std::sqrt (sum_squares);
#endif
			}
	};
	template<class T>
	struct matrix_vector_norm_inf:
		public matrix_vector_real_unary_functor<T> {
			typedef typename matrix_vector_real_unary_functor<T>::size_type size_type;
			typedef typename matrix_vector_real_unary_functor<T>::difference_type difference_type;
			typedef typename matrix_vector_real_unary_functor<T>::value_type value_type;
			typedef typename matrix_vector_real_unary_functor<T>::real_type real_type;
			typedef typename matrix_vector_real_unary_functor<T>::result_type result_type;
			typedef typename matrix_vector_real_unary_functor<T>::const_iterator1 const_iterator1;
			template<class E>
			static BOOST_UBLAS_INLINE
				result_type apply (size_type size,const const_iterator1 &it) {
					real_type t = real_type ();
					const_iterator1 it_temp = it;
					//size_type size (e ().size ());
					for (size_type i = 0; i < size; ++ i) {
						real_type u (std::abs (*it_temp));
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
						real_type u (std::abs (*it));
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
						real_type u (std::abs (*it));
						if (u > t) 
							t = u;
						++ it;
					}
					return t; 
			}
	};
	template<class E>
	struct matrix_vector_max:
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
					//size_type size (e ().size ());
					for (size_type i = 0; i < size; ++ i) {
						real_type u = *it_temp;
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
					real_type t = -(std::numeric_limits<result_type>::max());
					while (-- size >= 0) {
						real_type u = *it;
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
					real_type t = -(std::numeric_limits<result_type>::max());
					while (it != it_end) {
						real_type u = *it;
						if (u > t) 
							t = u;
						++ it;
					}
					return t; 
			}
	};
	template<class E>
	struct matrix_vector_rms:
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
						t += std::pow(*it_temp,2), ++ it_temp;
					return (std::sqrt(t));
			}
			// Dense case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (difference_type size, I it) {
					real_type t = real_type (0);
					while (-- size >= 0) {
						real_type u = std::pow(*it,2);
						t += u;
						++ it;
					}
					return (std::sqrt(t));
			}
			// Sparse case
			template<class I>
			static BOOST_UBLAS_INLINE
				result_type apply (I it, const I &it_end) { 
					real_type t = real_type ();
					while (it != it_end) {
						real_type u = std::pow(*it,2);
						t += u;
						++ it;
					}
					return (std::sqrt(t)); 
			}
	};
}}}
#endif //FUNCTIONAL_EXT_H
