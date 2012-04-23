#ifndef EXPRESSION_TYPES_EXT_H
#define EXPRESSION_TYPES_EXT_H

#include "stdafx.h"
#include "ublas_ext_config.hpp"
#include <boost/numeric/ublas/expression_types.hpp>
#include "matlab_tags.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#define GEN_TEML_SPECIALIZE_2(name)																\
	template<class U, class A>																	\
class vector_expression<name<U, A> >:														\
	public vector_expression_ext<name<U, A>, vect_ext_traits<name<U, A> > > {				\
public:																						\
	using vector_expression_ext<name<U, A>, vect_ext_traits<name<U, A> > >::operator ();	\
	using vector_expression_ext<name<U, A>, vect_ext_traits<name<U, A> > >::operator [];	\
	}

#define GEN_TEML_SPECIALIZE_1(name)														\
	template<class U>																	\
class vector_expression<name<U> >:													\
	public vector_expression_ext<name<U>, vect_ext_traits<name<U> > > {				\
public:																				\
	using vector_expression_ext<name<U>, vect_ext_traits<name<U> > >::operator ();	\
	using vector_expression_ext<name<U>, vect_ext_traits<name<U> > >::operator [];	\
	}

#define GEN_TEML_SPECIALIZE_MAT_1(name)													\
	template<class U>																	\
class matrix_expression<name<U> >:													\
	public matrix_expression_ext<name<U>, mat_ext_traits<name<U> > > {				\
public:																				\
	using matrix_expression_ext<name<U>, mat_ext_traits<name<U> > >::operator ();	\
	using matrix_expression_ext<name<U>, mat_ext_traits<name<U> > >::operator [];	\
	}

#define GEN_TEML_SPECIALIZE_MAT_2(name)													\
	template<class U, class A>															\
class matrix_expression<name<U, A> >:												\
	public matrix_expression_ext<name<U, A>, mat_ext_traits<name<U, A> > > {		\
public:																				\
	using matrix_expression_ext<name<U, A>, mat_ext_traits<name<U, A> > >::operator ();	\
	using matrix_expression_ext<name<U, A>, mat_ext_traits<name<U, A> > >::operator [];	\
	}

#define GEN_TEML_SPECIALIZE_MAT_3(name)													\
	template<class U, class A, class F>													\
class matrix_expression<name<U, A, F> >:											\
	public matrix_expression_ext<name<U, A, F>, mat_ext_traits<name<U, A, F> > > {	\
public:																				\
	using matrix_expression_ext<name<U, A, F>, mat_ext_traits<name<U, A, F> > >::operator ();	\
	using matrix_expression_ext<name<U, A, F>, mat_ext_traits<name<U, A, F> > >::operator [];	\
	}

namespace boost { namespace numeric { namespace ublas {
	//namespace  BNUE = boost::numeric::ublas_ext;
	template <class slice_type>
	BOOST_UBLAS_INLINE
		slice_type make_slice (typename slice_type::size_type start, typename slice_type::difference_type stride, typename slice_type::size_type stop) {
			BOOST_UBLAS_CHECK (stride, bad_index ());
			typename slice_type::difference_type  diff(stop - start);
			return slice_type (start, stride, (diff + ((stride > 0)? 1: -1))/stride);
	}

	// Simple Projections
	template<class V>
	BOOST_UBLAS_INLINE
		vector_slice<V> vect_subslice_ext (V &data, typename V::size_type start, typename V::difference_type stride, typename V::size_type stop) {
			typedef basic_slice<typename V::size_type, typename V::difference_type> slice_type;
			return vector_slice<V> (data, make_slice<slice_type> (start, stride, stop));
	}
	template<class V>
	BOOST_UBLAS_INLINE
		vector_slice<const V> vect_subslice_ext (const V &data, typename V::size_type start, typename V::difference_type stride, typename V::size_type stop)  {
			typedef basic_slice<typename V::size_type, typename V::difference_type> slice_type;
			return vector_slice<const V> (data, make_slice<slice_type> (start, stride, stop));
	}

	//Vector expression extentions
	template<class E, class V>
	// Vector extention traits
	class vector_expression_ext:
		public ublas_expression<E> {
	public:
		static const unsigned complexity = 0;
		typedef E expression_type;
		typedef vector_tag type_category;
		/* E can be an incomplete type - to define the following we would need more template arguments
		typedef typename E::size_type size_type;
		*/

		BOOST_UBLAS_INLINE
			const expression_type &operator () () const {
				return *static_cast<const expression_type *> (this);
		}
		BOOST_UBLAS_INLINE
			expression_type &operator () () {
				return *static_cast<expression_type *> (this);
		}


#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	private:
		// projection types
		typedef vector_range<E> vector_range_type;
		typedef vector_range<const E> const_vector_range_type;
		typedef vector_slice<E> vector_slice_type;
		typedef vector_slice<const E> const_vector_slice_type;
		// vector_indirect_type will depend on the A template parameter 
		typedef basic_range<> default_range;    // required to avoid range/slice name confusion
		typedef basic_slice<> default_slice;
		typedef const V const_vector_type;
		typedef V vector_type;
		typedef typename V::size_type size_type;
		typedef typename V::difference_type difference_type;
		typedef typename V::value_type value_type;
		typedef typename V::const_reference const_reference;
		typedef typename V::reference reference;

		template <class D>
		BOOST_UBLAS_INLINE
			size_type to_size(const last_elem<D>& offset) const{
				return (size_type)((operator ()().size() +  offset.get_offset()) - 1);
		}
		BOOST_UBLAS_INLINE
			size_type to_size(size_type offset) const{
				return offset;
		}
		BOOST_UBLAS_INLINE
			size_type get_size() const{
				return (size_type)(operator ()().size());
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size(const basic_slice_ext<T1, T2>& slice) const{
				return (size_type)( slice.bIsStartOffset_?((operator ()().size() +  slice.start_) - 1):slice.start_);
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size(const basic_slice_ext<T1, T2>& slice, bool bIsEnd) const{
				return (size_type)( slice.bIsEndOffset_?((operator ()().size() +  slice.end_) - 1):slice.end_);
		}	
	public:
		BOOST_UBLAS_INLINE
			const_vector_range_type operator () (const default_range &r) const {
				return const_vector_range_type (operator ()(), r);
		}
		BOOST_UBLAS_INLINE
			vector_range_type operator () (const default_range &r) {
				return vector_range_type (operator () (), r);
		}
		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (const default_slice &s) const {
				return const_vector_slice_type (operator () (), s);
		}
		BOOST_UBLAS_INLINE
			vector_slice_type operator () (const default_slice &s) {
				return vector_slice_type (operator () (), s);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			const vector_indirect<const E, indirect_array<A> > operator () (const indirect_array<A> &ia) const {
				return vector_indirect<const E, indirect_array<A> >  (operator () (), ia);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			vector_indirect<E, indirect_array<A> > operator () (const indirect_array<A> &ia) {
				return vector_indirect<E, indirect_array<A> > (operator () (), ia);
		}

		BOOST_UBLAS_INLINE
			const_vector_range_type project (const default_range &r) const {
				return const_vector_range_type (operator () (), r);
		}
		BOOST_UBLAS_INLINE
			vector_range_type project (const default_range &r) {
				return vector_range_type (operator () (), r);
		}
		BOOST_UBLAS_INLINE
			const_vector_slice_type project (const default_slice &s) const {
				return const_vector_slice_type (operator () (), s);
		}
		BOOST_UBLAS_INLINE
			vector_slice_type project (const default_slice &s) {
				return vector_slice_type (operator () (), s);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			const vector_indirect<const E, indirect_array<A> > project (const indirect_array<A> &ia) const {
				return vector_indirect<const E, indirect_array<A> > (operator () (), ia);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			vector_indirect<E, indirect_array<A> > project (const indirect_array<A> &ia) {
				return vector_indirect<E, indirect_array<A> > (operator () (), ia);
		}
		// Element access
		template <class D>
		BOOST_UBLAS_INLINE
			const_reference operator ()(const last_elem<D>& offset) const{
				return operator ()()(to_size(offset));
		}
		template <class D>
		BOOST_UBLAS_INLINE
			reference operator ()(const last_elem<D>& offset) {
				return operator ()()(to_size(offset));
		}
		template <class D>
		BOOST_UBLAS_INLINE
			const_reference operator [](const last_elem<D>& offset) const{
				return operator ()()(to_size(offset));
		}
		template <class D>
		BOOST_UBLAS_INLINE
			reference operator [](const last_elem<D>& offset) {
				return operator ()()(to_size(offset));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			vector_slice_type operator () (T1 start, T2 stop) const {
				return operator ()(start, 1, stop);
		}

		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (T1 start, T2 stop) {
				return operator ()(start, 1, stop);
		}

		BOOST_UBLAS_INLINE
			vector_slice_type operator () (const all_elem&) {
				return operator () (0, 1, get_size() - 1);
		}

		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (const all_elem&) const{
				return operator () (0, 1, get_size() - 1);
		}

		BOOST_UBLAS_INLINE
			vector_slice_type operator () (const rev_all&) {
				return operator () (get_size() - 1, -1, 0);
		}

		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (const rev_all&) const{
				return operator ()(get_size() - 1, -1, 0);
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (T1 start, difference_type stride, T2 stop) const {
				return vect_subslice_ext(operator ()(), to_size(start), stride, to_size(stop));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			vector_slice_type operator () (T1 start, difference_type stride, T2 stop)  {
				return vect_subslice_ext(operator ()(), to_size(start), stride, to_size(stop));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_vector_slice_type operator () (const basic_slice_ext<T1, T2>& slice) const {
				return vect_subslice_ext(operator ()(), to_size(slice), slice.stride_, to_size(slice));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			vector_slice_type operator () (const basic_slice_ext<T1, T2>& slice)  {
				return vect_subslice_ext(operator ()(), to_size(slice), slice.stride_, to_size(slice, true));
		}


#endif
	};



	template<class C>
	class vect_ext_traits {
	};

	template<class T, class A, template <class T1, class T2> class V>
	class vect_ext_traits<V<T, A> > {
	public:
		typedef typename A::size_type size_type;
		typedef typename A::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};

	template<class T>
	class vect_ext_traits<matrix_row<T> > {
	public:
		typedef typename T::size_type size_type;
		typedef typename T::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};

	template<class T>
	class vect_ext_traits<matrix_column<T> > {
	public:
		typedef typename T::size_type size_type;
		typedef typename T::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};

	template<class T>
	class vect_ext_traits<matrix_vector_slice<T> > {
	public:
		typedef typename T::size_type size_type;
		typedef typename T::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};

	template<class T>
	class vect_ext_traits<matrix_vector_range<T> > {
	public:
		typedef typename T::size_type size_type;
		typedef typename T::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};
	template<class E>
	class vect_ext_traits<vector_reference<E> > {
	public:
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename E::value_type value_type;
		typedef typename E::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<E>,
			typename E::const_reference,
			typename E::reference>::type reference;
	};

	template<class V>
	class vect_ext_traits<vector_slice<V> > {
	public:
		typedef const V const_vector_type;
		typedef V vector_type;
		typedef typename V::size_type size_type;
		typedef typename V::difference_type difference_type;
		typedef typename V::value_type value_type;
		typedef typename V::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<V>,
			typename V::const_reference,
			typename V::reference>::type reference;
	};

	template<class V>
	class vect_ext_traits<vector_range<V> > {
	public:
		typedef const V const_vector_type;
		typedef V vector_type;
		typedef typename V::size_type size_type;
		typedef typename V::difference_type difference_type;
		typedef typename V::value_type value_type;
		typedef typename V::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<V>,
			typename V::const_reference,
			typename V::reference>::type reference;
	};

	template<class E, class F>
	class vect_ext_traits<vector_unary<E, F> > {
	public:
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef typename boost::mpl::if_<boost::is_same<F, scalar_identity<value_type> >,
			typename E::reference,
			value_type>::type reference;
	};

	template<class U, class IA>
	class vect_ext_traits<vector_indirect<U, IA> > {
	public:
		typedef typename U::size_type size_type;
		typedef typename U::difference_type difference_type;
		typedef typename U::value_type value_type;
		typedef typename U::const_reference const_reference;
		typedef typename U::reference reference;
	};

	template<class E1, class E2, class F>
	class vect_ext_traits<vector_binary<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};

	template<class E1, class E2, class F>
	class vect_ext_traits<vector_binary_scalar1<E1, E2, F> > {
	public:
		typedef typename E2::size_type size_type;
		typedef typename E2::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};

	template<class E1, class E2, class F>
	class vect_ext_traits<vector_binary_scalar2<E1, E2, F> > {
	public:
		typedef typename E1::size_type size_type;
		typedef typename E1::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;

	};

	// Vector expression extention through specialization
	template<class E1, class E2, class F, template <class T1, class T2, class T3> class B>
	class vector_expression<B<E1, E2, F> >:
		public vector_expression_ext<B<E1, E2, F>, vect_ext_traits<B<E1, E2, F> > > {
	public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
		using vector_expression_ext<B<E1, E2, F>, vect_ext_traits<B<E1, E2, F> > >::operator ();
		using vector_expression_ext<B<E1, E2, F>, vect_ext_traits<B<E1, E2, F> > >::operator [];
#endif
	};


	GEN_TEML_SPECIALIZE_2(vector);	
	GEN_TEML_SPECIALIZE_2(vector_indirect);
	GEN_TEML_SPECIALIZE_2(vector_unary);
	GEN_TEML_SPECIALIZE_2(matrix_vector_indirect);
	GEN_TEML_SPECIALIZE_1(vector_reference);
	GEN_TEML_SPECIALIZE_1(vector_range);
	GEN_TEML_SPECIALIZE_1(vector_slice);
	GEN_TEML_SPECIALIZE_1(matrix_row);
	GEN_TEML_SPECIALIZE_1(matrix_column);
	GEN_TEML_SPECIALIZE_1(matrix_vector_slice);
	GEN_TEML_SPECIALIZE_1(matrix_vector_range);

	
	template<class T, class A>
	class vector_container<vector<T, A> >:
		public vector_expression<vector<T, A> > {
	public:
		static const unsigned complexity = 0;
		typedef vector<T, A>  container_type;
		typedef vector_tag type_category;

		BOOST_UBLAS_INLINE
			const container_type &operator () () const {
				return *static_cast<const container_type *> (this);
		}
		BOOST_UBLAS_INLINE
			container_type &operator () () {
				return *static_cast<container_type *> (this);
		}

#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
		using vector_expression<vector<T, A> >::operator ();
		using vector_expression<vector<T, A> >::operator [];
#endif
	};

	// Simple Projections
	template<class M>
	BOOST_UBLAS_INLINE
		matrix_slice<M> mat_subslice_ext (M &data, typename M::size_type start1, typename M::difference_type stride1, typename M::size_type size1, typename M::size_type start2, typename M::difference_type stride2, typename M::size_type size2) {
			typedef basic_slice<typename M::size_type, typename M::difference_type> slice_type;
			return matrix_slice<M> (data, make_slice<slice_type>  (start1, stride1, size1), make_slice<slice_type>  (start2, stride2, size2));
	}
	template<class M>
	BOOST_UBLAS_INLINE
		matrix_slice<const M> mat_subslice_ext (const M &data, typename M::size_type start1, typename M::difference_type stride1, typename M::size_type size1, typename M::size_type start2, typename M::difference_type stride2, typename M::size_type size2) {
			typedef basic_slice<typename M::size_type, typename M::difference_type> slice_type;
			return matrix_slice<const M> (data, make_slice<slice_type>  (start1, stride1, size1), make_slice<slice_type>  (start2, stride2, size2));
	}

	//Matrix expression extention
	template<class E, class M>
	class matrix_expression_ext:
		public ublas_expression<E> {
	public:
		static const unsigned complexity = 0;
		typedef E expression_type;
		typedef matrix_tag type_category;
		/* E can be an incomplete type - to define the following we would need more template arguments
		typedef typename E::size_type size_type;
		*/

		BOOST_UBLAS_INLINE
			const expression_type &operator () () const {
				return *static_cast<const expression_type *> (this);
		}
		BOOST_UBLAS_INLINE
			expression_type &operator () () {
				return *static_cast<expression_type *> (this);
		}

#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	private:
		// projection types
		typedef vector_range<E> vector_range_type;
		typedef const vector_range<const E> const_vector_range_type;
		typedef vector_slice<E> vector_slice_type;
		typedef const vector_slice<const E> const_vector_slice_type;
		typedef matrix_row<E> matrix_row_type;
		typedef const matrix_row<const E> const_matrix_row_type;
		typedef matrix_column<E> matrix_column_type;
		typedef const  matrix_column<const E> const_matrix_column_type;
		typedef matrix_range<E> matrix_range_type;
		typedef const matrix_range<const E> const_matrix_range_type;
		typedef matrix_slice<E> matrix_slice_type;
		typedef const matrix_slice<const E> const_matrix_slice_type;

		typedef vector_range<matrix_row< E > > matrix_row_range_type;
		typedef vector_range<matrix_row< const E > > const_matrix_row_range_type;
		typedef vector_range<matrix_column< E > > matrix_column_range_type;
		typedef vector_range<matrix_column< const E > > const_matrix_column_range_type;

		typedef vector_slice<matrix_row< E > >matrix_row_slice_type;
		typedef vector_slice<matrix_row< const E > > const_matrix_row_slice_type;
		typedef vector_slice<matrix_column< E > > matrix_column_slice_type;
		typedef vector_slice<matrix_column< const E > > const_matrix_column_slice_type;

		// matrix_indirect_type will depend on the A template parameter 
		typedef basic_range<> default_range;    // required to avoid range/slice name confusion
		typedef basic_slice<> default_slice;

		typedef const M const_matrix_type;
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::difference_type difference_type;
		typedef typename M::value_type value_type;
		typedef typename M::const_reference const_reference;
		typedef typename M::reference reference;

		template <class D>
		BOOST_UBLAS_INLINE
			size_type to_size1(const last_elem<D>& offset) const{
				return (size_type)((operator ()().size1() +  offset.get_offset()) - 1);
		}
		BOOST_UBLAS_INLINE
			size_type to_size1(size_type offset) const{
				return offset;
		}
		template <class D>
		BOOST_UBLAS_INLINE
			size_type to_size2(const last_elem<D>& offset) const{
				return (size_type)((operator ()().size2() +  offset.get_offset()) - 1);
		}
		BOOST_UBLAS_INLINE
			size_type to_size2(size_type offset) const{
				return offset;
		}		
		BOOST_UBLAS_INLINE
			size_type get_size1() const{
				return (size_type)(operator ()().size1());
		}		BOOST_UBLAS_INLINE
			size_type get_size2() const{
				return (size_type)(operator ()().size2());
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size1(const basic_slice_ext<T1, T2>& slice) const{
				return (size_type)( slice.bIsStartOffset_?((operator ()().size1() +  slice.start_) - 1):slice.start_);
		}		template <class T1, class T2>
			BOOST_UBLAS_INLINE
			size_type to_size2(const basic_slice_ext<T1, T2>& slice) const{
				return (size_type)( slice.bIsStartOffset_?((operator ()().size2() +  slice.start_) - 1):slice.start_);
		}		template <class T1, class T2>
			BOOST_UBLAS_INLINE
			size_type to_size1(const basic_slice_ext<T1, T2>& slice, bool bIsEnd) const{
				return (size_type)( slice.bIsEndOffset_?((operator ()().size1() +  slice.end_) - 1):slice.end_);
		}	
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size2(const basic_slice_ext<T1, T2>& slice, bool bIsEnd) const{
				return (size_type)(slice.bIsEndOffset_?((operator ()().size2() +  slice.end_) - 1):slice.end_);
		}	

		// For vector (linear) indexing of matrices
		template <class D>
		BOOST_UBLAS_INLINE
			size_type to_size(const last_elem<D>& offset) const{
				return (size_type)((get_size() +  offset.get_offset()) - 1);
		}
		BOOST_UBLAS_INLINE
			size_type to_size(size_type offset) const{
				return offset;
		}
		BOOST_UBLAS_INLINE
			size_type get_size() const{
				return (size_type)((operator ()().size1())*(operator ()().size2()));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size(const basic_slice_ext<T1, T2>& slice) const{
				return (size_type)( slice.bIsStartOffset_?((get_size() +  slice.start_) - 1):slice.start_);
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			size_type to_size(const basic_slice_ext<T1, T2>& slice, bool bIsEnd) const{
				return (size_type)( slice.bIsEndOffset_?((get_size() +  slice.end_) - 1):slice.end_);
		}		public:
		// Element access
		template <class D>
		BOOST_UBLAS_INLINE
			const_reference operator ()(const last_elem<D>& offset1,const last_elem<D>& offset2 ) const{
				return operator ()()(to_size1(offset1), to_size2(offset2));
		}
		template <class D>
		BOOST_UBLAS_INLINE
			reference operator ()(const last_elem<D>& offset1,const last_elem<D>& offset2 ) {
				return operator ()()(to_size1(offset1), to_size2(offset2));
		}

		template <class T>
		BOOST_UBLAS_INLINE
			const_matrix_row_type operator [] (T offset) const {
				return const_matrix_row_type (operator () (), to_size1(offset));
		}
		template <class T>
		BOOST_UBLAS_INLINE
			matrix_row_type operator [] (T offset) {
				return matrix_row_type (operator () (), to_size1(offset));
		}
		template <class T>
		BOOST_UBLAS_INLINE
			const_matrix_row_type row (T offset) const {
				return const_matrix_row_type (operator () (),  to_size1(offset));
		}
		template <class T>
		BOOST_UBLAS_INLINE
			matrix_row_type row (T offset) {
				return matrix_row_type (operator () (),  to_size1(offset));
		}
		template <class T>
		BOOST_UBLAS_INLINE
			const_matrix_column_type column (T offset) const {
				return const_matrix_column_type (operator () (),  to_size1(offset));
		}
		template <class T>
		BOOST_UBLAS_INLINE
			matrix_column_type column (T offset) {
				return matrix_column_type (operator () (),  to_size1(offset));
		}

		BOOST_UBLAS_INLINE
			const_matrix_range_type operator () (const default_range &r1, const default_range &r2) const {
				return const_matrix_range_type (operator () (), r1, r2);
		}
		BOOST_UBLAS_INLINE
			matrix_range_type operator () (const default_range &r1, const default_range &r2) {
				return matrix_range_type (operator () (), r1, r2);
		}
		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const default_slice &s1, const default_slice &s2) const {
				return const_matrix_slice_type (operator () (), s1, s2);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const default_slice &s1, const default_slice &s2) {
				return matrix_slice_type (operator () (), s1, s2);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			const matrix_indirect<const E, indirect_array<A> > operator () (const indirect_array<A> &ia1, const indirect_array<A> &ia2) const {
				return matrix_indirect<const E, indirect_array<A> > (operator () (), ia1, ia2);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			matrix_indirect<E, indirect_array<A> > operator () (const indirect_array<A> &ia1, const indirect_array<A> &ia2) {
				return matrix_indirect<E, indirect_array<A> > (operator () (), ia1, ia2);
		}

		BOOST_UBLAS_INLINE
			const_matrix_range_type project (const default_range &r1, const default_range &r2) const {
				return const_matrix_range_type (operator () (), r1, r2);
		}
		BOOST_UBLAS_INLINE
			matrix_range_type project (const default_range &r1, const default_range &r2) {
				return matrix_range_type (operator () (), r1, r2);
		}
		BOOST_UBLAS_INLINE
			const_matrix_slice_type project (const default_slice &s1, const default_slice &s2) const {
				return const_matrix_slice_type (operator () (), s1, s2);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type project (const default_slice &s1, const default_slice &s2) {
				return matrix_slice_type (operator () (), s1, s2);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			const matrix_indirect<const E, indirect_array<A> > project (const indirect_array<A> &ia1, const indirect_array<A> &ia2) const {
				return matrix_indirect<const E, indirect_array<A> > (operator () (), ia1, ia2);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			matrix_indirect<E, indirect_array<A> > project (const indirect_array<A> &ia1, const indirect_array<A> &ia2) {
				return matrix_indirect<E, indirect_array<A> > (operator () (), ia1, ia2);
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const basic_slice_ext<T1, T2>& slice1, const basic_slice_ext<T1, T2>& slice2) const {
				return mat_subslice_ext(operator ()(), to_size1(slice1), slice1.stride_, to_size1(slice1, true), to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const basic_slice_ext<T1, T2>& slice1, const basic_slice_ext<T1, T2>& slice2) {
				return mat_subslice_ext(operator ()(), to_size1(slice1), slice1.stride_, to_size1(slice1, true), to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}

		template <class T1, class T2, class T3>
		BOOST_UBLAS_INLINE
			const_matrix_row_slice_type operator () (const last_elem<T3>& offset1, const basic_slice_ext<T1, T2>& slice2) const {
				return row(to_size1(offset1))(to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}
		template <class T1, class T2, class T3>
		BOOST_UBLAS_INLINE
			matrix_row_slice_type operator () (const last_elem<T3>& offset1, const basic_slice_ext<T1, T2>& slice2) {
				return row(to_size1(offset1))(to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_matrix_row_slice_type operator () (size_type offset1, const basic_slice_ext<T1, T2>& slice2) const {
				return row(to_size1(offset1))(to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			matrix_row_slice_type operator () (size_type offset1, const basic_slice_ext<T1, T2>& slice2) {
				return row(to_size1(offset1))(to_size2(slice2), slice2.stride_, to_size2(slice2, true));
		}


		template <class T3>
		BOOST_UBLAS_INLINE
			const_matrix_row_slice_type operator () (const last_elem<T3>& offset1, const all_elem&) const {
				return row(to_size1(offset1))(0, 1, get_size2() - 1);
		}
		template <class T3>
		BOOST_UBLAS_INLINE
			matrix_row_slice_type operator () (const last_elem<T3>& offset1, const all_elem&) {
				return row(to_size1(offset1))(0, 1, get_size2() - 1);
		}

		template <class T3>
		BOOST_UBLAS_INLINE
			const_matrix_row_slice_type operator () (const last_elem<T3>& offset1, const rev_all&) const {
				return row(to_size1(offset1))(get_size2() - 1, -1, 0);
		}
		template <class T3>
		BOOST_UBLAS_INLINE
			matrix_row_slice_type operator () (const last_elem<T3>& offset1, const rev_all&) {
				return row(to_size1(offset1))(get_size2() - 1, -1, 0);
		}

		template <class T1, class T2, class T3>
		BOOST_UBLAS_INLINE
			const_matrix_column_slice_type operator () (const basic_slice_ext<T1, T2>& slice1, const last_elem<T3>& offset2) const {
				return column(to_size1(offset2))( to_size1(slice1), slice1.stride_, to_size1(slice1, true));
		}
		template <class T1, class T2, class T3>
		BOOST_UBLAS_INLINE
			matrix_column_slice_type operator ()(const basic_slice_ext<T1, T2>& slice1, const last_elem<T3>& offset2) {
				return column(to_size1(offset2))(to_size1(slice1), slice1.stride_, to_size1(slice1, true));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			const_matrix_column_slice_type operator () (const basic_slice_ext<T1, T2>& slice1, size_type offset2) const {
				return column(to_size1(offset2))( to_size1(slice1), slice1.stride_, to_size1(slice1, true));
		}
		template <class T1, class T2>
		BOOST_UBLAS_INLINE
			matrix_column_slice_type operator ()(const basic_slice_ext<T1, T2>& slice1, size_type offset2) {
				return column(to_size1(offset2))(to_size1(slice1), slice1.stride_, to_size1(slice1, true));
		}

		template <class T3>
		BOOST_UBLAS_INLINE
			const_matrix_column_slice_type operator () (const all_elem&, const last_elem<T3>& offset2) const {
				return column(to_size1(offset2))(0, 1, get_size1() - 1);
		}
		template <class T3>
		BOOST_UBLAS_INLINE
			matrix_column_slice_type operator ()(const all_elem&, const last_elem<T3>& offset2) {
				return column(to_size1(offset2))(0, 1, get_size1() - 1);
		}

		template <class T3>
		BOOST_UBLAS_INLINE
			const_matrix_column_slice_type operator () (const rev_all&, const last_elem<T3>& offset2) const {
				return column(to_size1(offset2))(get_size1() - 1, -1, 0);
		}
		template <class T3>
		BOOST_UBLAS_INLINE
			matrix_column_slice_type operator ()(const rev_all&, const last_elem<T3>& offset2) {
				return column(to_size1(offset2))(get_size1() - 1, -1, 0);
		}

		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const all_elem&, const all_elem&) const {
				return mat_subslice_ext(operator ()(), 0, 1, get_size1() - 1, 0, 1, get_size2() - 1);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const all_elem&, const all_elem&) {
				return mat_subslice_ext(operator ()(), 0, 1, get_size1() - 1, 0, 1, get_size2() - 1);
		}

		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const rev_all&, const rev_all&) const {
				return mat_subslice_ext(operator ()(), get_size1() - 1, -1, 0, get_size2() - 1, -1, 0);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const rev_all&, const rev_all&) {
				return mat_subslice_ext(operator ()(), get_size1() - 1, -1, 0, get_size2() - 1, -1, 0);
		}

		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const all_elem&, const rev_all&) const {
				return mat_subslice_ext(operator ()(), 0, 1, get_size1() - 1,get_size2() - 1, -1, 0);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const all_elem&, const rev_all&) {
				return mat_subslice_ext(operator ()(), 0, 1, get_size1() - 1, get_size2() - 1, -1, 0);
		}

		BOOST_UBLAS_INLINE
			const_matrix_slice_type operator () (const rev_all&, const all_elem&) const {
				return mat_subslice_ext(operator ()(), get_size1() - 1, -1, 0, 0, 1, get_size2() - 1);
		}
		BOOST_UBLAS_INLINE
			matrix_slice_type operator () (const rev_all&, const all_elem&) {
				return mat_subslice_ext(operator ()(), get_size1() - 1, -1, 0, 0, 1, get_size2() - 1);
		}
		template<class A, class T>
		BOOST_UBLAS_INLINE
			const vector_indirect<matrix_row<const E>, indirect_array<A> >  operator () (const last_elem<T>& offset1, const indirect_array<A> &ia2) const {
				return operator()(to_size1(offset1), ia2);

		}
		template<class A, class T>
		BOOST_UBLAS_INLINE
			vector_indirect< matrix_row<E>, indirect_array<A> >  operator () (const last_elem<T>& offset1, const indirect_array<A> &ia2)  {
				return operator()(to_size1(offset1), ia2);
		}

		template<class A>
		BOOST_UBLAS_INLINE
			const vector_indirect<matrix_row<const E>, indirect_array<A> >  operator () (size_type offset1, const indirect_array<A> &ia2) const {
				return vector_indirect<matrix_row<const E>, indirect_array<A> > (row(offset1), ia2);

		}
		template<class A>
		BOOST_UBLAS_INLINE
			vector_indirect< matrix_row<E>, indirect_array<A> >  operator () (size_type offset1, const indirect_array<A> &ia2)  {
				return vector_indirect<matrix_row<E>, indirect_array<A> > (row(offset1), ia2);
		}

		template<class A, class T>
		BOOST_UBLAS_INLINE
			const vector_indirect<matrix_column<const E>, indirect_array<A> >  operator () (const indirect_array<A> &ia1, const last_elem<T>& offset2) const {
				return operator () (ia1, to_size1(offset2));
		}
		template<class A, class T>
		BOOST_UBLAS_INLINE
			vector_indirect<matrix_column<E>, indirect_array<A> >  operator () (const indirect_array<A> &ia1, const last_elem<T>& offset2) {
				return operator () (ia1, to_size1(offset2));
		}

		template<class A>
		BOOST_UBLAS_INLINE
			const vector_indirect<matrix_column<const E>, indirect_array<A> >  operator () (const indirect_array<A> &ia1, size_type offset2) const {
				return vector_indirect<matrix_column<const E>, indirect_array<A> > (column(offset2), ia1);
		}
		template<class A>
		BOOST_UBLAS_INLINE
			vector_indirect<matrix_column<E>, indirect_array<A> >  operator () (const indirect_array<A> &ia1, size_type offset2) {
				return vector_indirect<matrix_column<E>, indirect_array<A> > (column(offset2), ia1);
		}


		// Vector indexing of matrix TODO
		//template<class A>
		//BOOST_UBLAS_INLINE
		//matrix_vector_indirect<const E, indirect_array<A> >  operator () (const indirect_array<A>& ia) const {
		//		return matrix_vector_indirect<const E, indirect_array<A> >(operator ()(), ia);
		//}
		//template<class A>
		//BOOST_UBLAS_INLINE
		//matrix_vector_indirect<E, indirect_array<A> >  operator () (const indirect_array<A>& ia)  {
		//		return matrix_vector_indirect<const E, indirect_array<A> >(operator ()(), ia);
		//}

#endif
	};


	// Matrix extention traits
	template<class C>
	class mat_ext_traits {
	};

	template<class T, class L, class A, template <class T1, class T2, class T3> class M>
	class mat_ext_traits<M<T, L, A> > {
	public:
		typedef typename A::size_type size_type;
		typedef typename A::difference_type difference_type;
		typedef T value_type;
		typedef typename type_traits<T>::const_reference const_reference;
		typedef T &reference;
	};

	template<class E>	class mat_ext_traits<matrix_reference<E> > {
	public:
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename E::value_type value_type;
		typedef typename E::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<E>,
			typename E::const_reference,
			typename E::reference>::type reference;
	};

	template<class V>
	class mat_ext_traits<matrix_slice<V> > {
	public:
		typedef const V const_matrix_type;
		typedef V matrix_type;
		typedef typename V::size_type size_type;
		typedef typename V::difference_type difference_type;
		typedef typename V::value_type value_type;
		typedef typename V::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<V>,
			typename V::const_reference,
			typename V::reference>::type reference;
	};

	template<class V>
	class mat_ext_traits<matrix_range<V> > {
	public:
		typedef const V const_matrix_type;
		typedef V matrix_type;
		typedef typename V::size_type size_type;
		typedef typename V::difference_type difference_type;
		typedef typename V::value_type value_type;
		typedef typename V::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<V>,
			typename V::const_reference,
			typename V::reference>::type reference;
	};

	template<class E, class F>
	class mat_ext_traits<matrix_unary1<E, F> > {
	public:
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};
	template<class E, class F>
	class mat_ext_traits<matrix_unary2<E, F> > {
	public:
		typedef typename E::size_type size_type;
		typedef typename E::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef typename boost::mpl::if_<boost::is_same<F, scalar_identity<value_type> >,
			typename E::reference,
			value_type>::type reference;
	};


	template<class U, class IA>
	class mat_ext_traits<matrix_indirect<U, IA> > {
	public:
		typedef typename U::size_type size_type;
		typedef typename U::difference_type difference_type;
		typedef typename U::value_type value_type;
		typedef typename U::const_reference const_reference;
		typedef typename U::reference reference;
	};

	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_binary<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};

	template<class E1, class E2, class F>
	class mat_ext_traits<vector_matrix_binary<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};
	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_binary_scalar1<E1, E2, F> > {
	public:
		typedef typename E2::size_type size_type;
		typedef typename E2::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;
	};
	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_binary_scalar2<E1, E2, F> > {
	public:
		typedef typename E1::size_type size_type;
		typedef typename E1::difference_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;

	};
	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_vector_binary1<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;

	};
	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_vector_binary2<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;

	};
	template<class E1, class E2, class F>
	class mat_ext_traits<matrix_matrix_binary<E1, E2, F> > {
	public:
		typedef typename promote_traits<typename E1::size_type, typename E2::size_type>::promote_type size_type;
		typedef typename promote_traits<typename E1::difference_type, typename E2::difference_type>::promote_type difference_type;
		typedef typename F::result_type value_type;
		typedef value_type const_reference;
		typedef const_reference reference;

	};



	// Matrix expression extention through specialization

	GEN_TEML_SPECIALIZE_MAT_1(matrix_range);
	GEN_TEML_SPECIALIZE_MAT_1(matrix_slice);
	GEN_TEML_SPECIALIZE_MAT_1(matrix_reference);
	GEN_TEML_SPECIALIZE_MAT_2(matrix_indirect);
	GEN_TEML_SPECIALIZE_MAT_2(matrix_unary1);
	GEN_TEML_SPECIALIZE_MAT_2(matrix_unary2);
	GEN_TEML_SPECIALIZE_MAT_3(matrix);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_binary);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_vector_binary1);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_vector_binary2);
	GEN_TEML_SPECIALIZE_MAT_3(vector_matrix_binary);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_binary_scalar1);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_binary_scalar2);
	GEN_TEML_SPECIALIZE_MAT_3(matrix_matrix_binary);

	template<class T, class L, class A>
	class matrix_container<matrix<T, L, A> >:
		public matrix_expression<matrix<T, L, A> > {
	public:
		static const unsigned complexity = 0;
		typedef matrix<T, L, A>  container_type;
		typedef matrix_tag type_category;

		BOOST_UBLAS_INLINE
			const container_type &operator () () const {
				return *static_cast<const container_type *> (this);
		}
		BOOST_UBLAS_INLINE
			container_type &operator () () {
				return *static_cast<container_type *> (this);
		}

#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
		using matrix_expression<matrix<T, L, A> >::operator ();
		using matrix_expression<matrix<T, L, A> >::operator [];
#endif
	};

	typedef indirect_array<> index_vector;

}}}

#endif

