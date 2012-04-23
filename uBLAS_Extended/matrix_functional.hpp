#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/functional.hpp>
#include <cmath>
#include <cstdlib>

using namespace boost::numeric::ublas;

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
                real_type u (type_traits<value_type>::norm_1 (*it_temp));
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
                real_type u (type_traits<value_type>::norm_1 (*it));
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
                real_type u (type_traits<value_type>::norm_1 (*it));
                t += u;
                ++ it;
            }
            return t;
        }
    };



//norm_L1 m[j] = norm_L1 (m[i][j])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_vector_unary_traits<E, matrix_vector_norm_1 <E> >::result_type
norm_L1 (const matrix_expression<E> &e) {
	  typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_1 <E> >::expression_type expression_type;
		  return expression_type (e ());
}




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
                real_type u (type_traits<value_type>::norm_2 (*it_temp));
                t +=  u * u, ++ it_temp;
						}
            return type_traits<real_type>::sqrt (t);
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            for (size_type i = 0; i < size; ++ i) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                if (scale < u) {
                    real_type v (scale / u);
                    sum_squares = sum_squares * v * v + real_type (1);
                    scale = u;
                } else {
                    real_type v (u / scale);
                    sum_squares += v * v;
                }
            }
            return scale * type_traits<real_type>::sqrt (sum_squares);
#endif
        }
        // Dense case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (difference_type size, I it) {
#ifndef BOOST_UBLAS_SCALED_NORM
            real_type t = real_type ();
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return type_traits<real_type>::sqrt (t);
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            while (-- size >= 0) {
                real_type u (type_traits<value_type>::norm_2 (*it));
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
            return scale * type_traits<real_type>::sqrt (sum_squares);
#endif
        }
        // Sparse case
        template<class I>
        static BOOST_UBLAS_INLINE
        result_type apply (I it, const I &it_end) {
#ifndef BOOST_UBLAS_SCALED_NORM
            real_type t = real_type ();
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_2 (*it));
                t +=  u * u;
                ++ it;
            }
            return type_traits<real_type>::sqrt (t);
#else
            real_type scale = real_type ();
            real_type sum_squares (1);
            while (it != it_end) {
                real_type u (type_traits<value_type>::norm_2 (*it));
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
            return scale * type_traits<real_type>::sqrt (sum_squares);
#endif
        }
    };




//norm_L2 m[j] = norm_L2 (m[i][j])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_vector_unary_traits<E, matrix_vector_norm_2 <E> >::result_type
norm_L2 (const matrix_expression<E> &e) {
	    typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_2 <E> >::expression_type expression_type;
			return expression_type (e ());
}




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
                real_type u (type_traits<value_type>::norm_inf (*it_temp));
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
                real_type u (type_traits<value_type>::norm_inf (*it));
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
                real_type u (type_traits<value_type>::norm_inf (*it));
                if (u > t) 
                    t = u;
                ++ it;
            }
            return t; 
        }
    };




//norm_inf m[j] = norm_inf (m[i][j])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_vector_unary_traits<E, matrix_vector_norm_inf <E> >::result_type
col_norm_inf (const matrix_expression<E> &e) {
	    typedef typename matrix_vector_unary_traits<E, matrix_vector_norm_inf <E> >::expression_type expression_type;
			 return expression_type (e ());
}




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


// (m ^ t) [i] [j] = m [i] [j] ^ t
template<class E1, class T2>
BOOST_UBLAS_INLINE
typename matrix_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::result_type
operator ^ (const matrix_expression<E1> &e1, const T2 &e2) {
	  typedef typename matrix_binary_scalar2_traits<E1, const T2, scalar_power<typename E1::value_type, T2> >::expression_type expression_type;
		 return expression_type (e1(), e2);
}


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


// (log10 v) [i] = log10 (v [i])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_log10<typename E::value_type> >::result_type
log10 (const matrix_expression<E> &e) {
      typedef typename matrix_unary1_traits<E, scalar_log10<typename E::value_type> >::expression_type expression_type;
      return expression_type (e ());
}



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


// (abs v) [i] = abs (v [i])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_abs<typename E::value_type> >::result_type
abs (const matrix_expression<E> &e) {
    typedef typename matrix_unary1_traits<E, scalar_abs<typename E::value_type> >::expression_type expression_type;
    return expression_type (e ());
}



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




// (exp v) [i] = exp (v [i])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::result_type
exp (const matrix_expression<E> &e) {
     typedef typename matrix_unary1_traits<E, scalar_exp<typename E::value_type> >::expression_type expression_type;
     return expression_type (e ());
}



template<class T>
struct scalar_log:
public scalar_unary_functor<T> {
  typedef typename scalar_unary_functor<T>::argument_type argument_type;
  typedef typename scalar_unary_functor<T>::result_type result_type;

  static BOOST_UBLAS_INLINE
  result_type apply (argument_type t) {
      return std::log(t);
  }
};



// (log v) [i] = log (v [i])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::result_type
log (const matrix_expression<E> &e) {
     typedef typename matrix_unary1_traits<E, scalar_log<typename E::value_type> >::expression_type expression_type;
     return expression_type (e ());
}

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




// (sqrt v) [i] = sqrt (v [i])
template<class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits<E, scalar_sqrt<typename E::value_type> >::result_type
sqrt  (const matrix_expression<E> &e) {
      typedef typename matrix_unary1_traits<E, scalar_sqrt<typename E::value_type> >::expression_type expression_type;
      return expression_type (e ());
}




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




// (sin v) [i] = sin (v [i])
 template<class E>
 BOOST_UBLAS_INLINE
 typename matrix_unary1_traits<E, scalar_sin<typename E::value_type> >::result_type
 sin  (const matrix_expression<E> &e) {
       typedef typename matrix_unary1_traits<E, scalar_sin<typename E::value_type> >::expression_type expression_type;
       return expression_type (e ());
 }





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




// (cos v) [i] = cos (v [i])
 template<class E>
 BOOST_UBLAS_INLINE
 typename matrix_unary1_traits<E, scalar_cos<typename E::value_type> >::result_type
 cos  (const matrix_expression<E> &e) {
       typedef typename matrix_unary1_traits<E, scalar_cos<typename E::value_type> >::expression_type expression_type;
       return expression_type (e ());
 }



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



template <class E>
BOOST_UBLAS_INLINE
typename matrix_unary1_traits <E, logical_not < typename E::value_type> > :: result_type
operator ! (const matrix_expression <E> &e) {
  typedef typename matrix_unary1_traits<E, logical_not <typename E::value_type> > :: expression_type expression_type;
  return expression_type (e () );
}



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



// ((m > t) [i][j]  =  m [i][j] > t
 template<class E1>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar2_traits<E1,typename E1::value_type, scalar_greater<typename E1::value_type, typename E1::value_type> >::result_type
 operator > (const matrix_expression<E1> &e1, typename E1::value_type e2) {
        typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
            return expression_type (e1 (), e2);
}



 // (t > m) [i][j]  = t > m [i][j]
 template<class E2>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater<typename E2::value_type, typename E2::value_type> >::result_type
 operator > (typename E2::value_type e1, const matrix_expression<E2> &e2) {
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



// (m >= t) [i][j] = m [i][j] >= t
 template<class E1>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::result_type
 operator >= (const matrix_expression<E1> &e1, typename E1::value_type e2) {
        typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_greater_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
            return expression_type (e1 (), e2);
}



 // (t >= m) [i][j] = t >= v [i][j]
 template<class E2>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_greater_equal<typename E2::value_type, typename E2::value_type> >::result_type
 operator >= (typename E2::value_type e1, const matrix_expression<E2> &e2) {
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



// (m < t) [i][j] = m [i][j] < t
 template<class E1>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser<typename E1::value_type, typename E1::value_type> >::result_type
 operator < (const matrix_expression<E1> &e1, typename E1::value_type e2) {
        typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
            return expression_type (e1 (), e2);
}



 // (m < v) [i][j] = t < m [i][j]
 template<class E2>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser<typename E2::value_type, typename E2::value_type> >::result_type
 operator < (typename E2::value_type e1, const matrix_expression<E2> &e2) {
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



// (m <= t) [i][j] = m [i][j] <= t
 template<class E1>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::result_type
 operator <= (const matrix_expression<E1> &e1, typename E1::value_type e2) {
        typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_lesser_equal<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
            return expression_type (e1 (), e2);
}



 // (t <= m) [i][j] = t <= m [i][j]
 template<class E2>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_lesser_equal<typename E2::value_type, typename E2::value_type> >::result_type
 operator <= (typename E2::value_type e1, const matrix_expression<E2> &e2) {
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



// (m != t) [i][j] = m [i][j] != t
 template<class E1>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_inequality<typename E1::value_type, typename E1::value_type> >::result_type
 operator != (const matrix_expression<E1> &e1, typename E1::value_type e2) {
        typedef typename matrix_binary_scalar2_traits<E1, typename E1::value_type, scalar_inequality<typename E1::value_type, typename E1::value_type> >::expression_type expression_type;
            return expression_type (e1 (), e2);
}



 // (t != m) [i][j] = t != m [i][j]
 template<class E2>
 BOOST_UBLAS_INLINE
 typename matrix_binary_scalar1_traits<typename E2::value_type, E2, scalar_inequality<typename E2::value_type, typename E2::value_type> >::result_type
 operator != (typename E2::value_type e1, const matrix_expression<E2> &e2) {
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

