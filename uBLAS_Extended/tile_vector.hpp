#ifndef _BOOST_TILE_MATRIX_FROM_VECTOR_
#define _BOOST_TILE_MATRIX_FROM_VECTOR_

// Created by Karthick Manivannan and Raghavender Boorgapally

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/detail/temporary.hpp>

namespace boost { namespace numeric { namespace ublas {


	template<class M>
	class vector_tile :
		public matrix_expression<vector_tile<M> > {

			typedef vector_tile<M> self_type;
	public:
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
		using matrix_expression<self_type>::operator ();
#endif
		typedef const M const_matrix_type;
		typedef M matrix_type;
		typedef typename M::size_type size_type;
		typedef typename M::difference_type difference_type;
		typedef typename M::value_type value_type;
		typedef typename M::const_reference const_reference;
		typedef typename boost::mpl::if_<boost::is_const<M>,
			typename M::const_reference,
			typename M::reference>::type reference;
		typedef typename boost::mpl::if_<boost::is_const<M>,
			typename M::const_closure_type,
			typename M::closure_type>::type matrix_closure_type;
		typedef const self_type const_closure_type;
		typedef self_type closure_type;
		// Replaced by _temporary_traits to avoid type requirements on M
		//typedef typename M::vector_temporary_type vector_temporary_type;
		//typedef typename M::matrix_temporary_type matrix_temporary_type;
		typedef typename storage_restrict_traits<typename M::storage_category,
			packed_proxy_tag>::storage_category storage_category;



		// Construction and destruction
		BOOST_UBLAS_INLINE
			vector_tile(matrix_type &data):
		matrix_expression<self_type> (),
			data_ (data), num_rows (1), num_cols (data.size()) {}

		BOOST_UBLAS_INLINE
			vector_tile (const vector_tile &m):
		matrix_expression<self_type> (),
			data_ (m.data_), num_rows (m.num_rows), num_cols (m.num_cols) {}

		BOOST_UBLAS_INLINE
			vector_tile (matrix_type &data, size_type rows, size_type cols):
		matrix_expression<self_type> (),
			data_ (data), num_rows (rows), num_cols (cols*data.size()) {}



		// Accessors
		BOOST_UBLAS_INLINE
			size_type size1 () const {
				return num_rows;
		}
		BOOST_UBLAS_INLINE
			size_type size2 () const {
				return num_cols;
		}


		// Storage accessors
		BOOST_UBLAS_INLINE
			const matrix_closure_type &data () const {
				return data_;
		}
		BOOST_UBLAS_INLINE
			matrix_closure_type &data () {
				return data_;
		}



		// Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
		BOOST_UBLAS_INLINE
			const_reference operator () (size_type i, size_type j) const {
				BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
				BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
				//std::cout << "(" << i << "," << j << ")" << std::endl;
				return data () (j%data_.size());
		}
		BOOST_UBLAS_INLINE
			reference operator () (size_type i, size_type j) {
				BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
				BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
				//std::cout << "(" << i << "," << j << ")" << std::endl;
				return data () (j%data_.size());
		}
#else
		BOOST_UBLAS_INLINE
			reference operator () (size_type i, size_type j) const {
				BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
				BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
				//std::cout << "(" << i << "," << j << ")" << std::endl;
				return data () (j%data_.size());
		}
#endif

		// Assignment
		// Not allowed!

        // Closure comparison
        BOOST_UBLAS_INLINE
        bool same_closure (const vector_tile &ba) const {
            return (*this).data ().same_closure (ba.data ());
        }
		// Iterator types
	private:
		// Use the matrix iterator
		typedef typename M::const_iterator const_subiterator1_type;
		typedef typename boost::mpl::if_<boost::is_const<M>,
			typename M::const_iterator,
			typename M::iterator>::type subiterator1_type;
		typedef typename M::const_iterator const_subiterator2_type;
		typedef typename boost::mpl::if_<boost::is_const<M>,
			typename M::const_iterator,
			typename M::iterator>::type subiterator2_type;

	public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
		typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
		typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
		typedef indexed_const_iterator1<self_type, packed_random_access_iterator_tag> const_iterator1;
		typedef indexed_const_iterator2<self_type, packed_random_access_iterator_tag> const_iterator2;
#else
		class const_iterator1;
		class iterator1;
		class const_iterator2;
		class iterator2;
#endif
		typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
		typedef reverse_iterator_base1<iterator1> reverse_iterator1;
		typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
		typedef reverse_iterator_base2<iterator2> reverse_iterator2;

		// Element lookup
		BOOST_UBLAS_INLINE
			const_iterator1 find1 (int rank, size_type i, size_type j) const {
				return const_iterator1 (*this, i ,j);
		}
		BOOST_UBLAS_INLINE
			iterator1 find1 (int rank, size_type i, size_type j) {
				return iterator1 (*this, i, j );
		}
		BOOST_UBLAS_INLINE
			const_iterator2 find2 (int rank, size_type i, size_type j) const {
				return const_iterator2 (*this, i, j);
		}
		BOOST_UBLAS_INLINE
			iterator2 find2 (int rank, size_type i, size_type j) {
				return iterator2 (*this, i, j);
		}

		// Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
		class const_iterator1:
			public container_const_reference<vector_tile>,
			public random_access_iterator_base<typename iterator_restrict_traits<
			typename const_subiterator1_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
			const_iterator1, value_type> {
		public:
			typedef typename const_subiterator1_type::value_type value_type;
			typedef typename const_subiterator1_type::difference_type difference_type;
			typedef typename const_subiterator1_type::reference reference;
			typedef typename const_subiterator1_type::pointer pointer;

			typedef const_iterator2 dual_iterator_type;
			typedef const_reverse_iterator2 dual_reverse_iterator_type;

			// Construction and destruction
			BOOST_UBLAS_INLINE
				const_iterator1 ():
			container_const_reference<self_type> (), rownum (), colnum () {}
			BOOST_UBLAS_INLINE
				const_iterator1 (const self_type &m, const size_type row_number, const size_type column_number):
			container_const_reference<self_type> (m), rownum (row_number), colnum (column_number) {}
			BOOST_UBLAS_INLINE
				const_iterator1 (const iterator1 &it):
			container_const_reference<self_type> (it ()), rownum(it.rownum), colnum(it.colnum) {}

			// Arithmetic
			BOOST_UBLAS_INLINE
				const_iterator1 &operator ++ () {
					++ rownum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator1 &operator -- () {
					--rownum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator1 &operator += (difference_type n) {
					rownum += n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator1 &operator -= (difference_type n) {
					rownum -= n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				difference_type operator - (const const_iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (rownum - it.rownum);
			}

			// Dereference
			BOOST_UBLAS_INLINE
				const_reference operator * () const {
					size_type i = index1 ();
					size_type j = index2 ();
					BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
					BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
					return (*this) () (i, j);
			}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_iterator2 begin () const {
					return (*this) ().find2 (1, index1 (), 0);
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_iterator2 end () const {
					return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_reverse_iterator2 rbegin () const {
					return const_reverse_iterator2 (end ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_reverse_iterator2 rend () const {
					return const_reverse_iterator2 (begin ());
			}
#endif

			// Indices
			BOOST_UBLAS_INLINE
				size_type index1 () const {
					return rownum;
			}
			BOOST_UBLAS_INLINE
				size_type index2 () const {
					return colnum;
			}

			// Assignment
			BOOST_UBLAS_INLINE
				const_iterator1 &operator = (const const_iterator1 &it) {
					container_const_reference<self_type>::assign (&it ());
					rownum = it.rownum;
					colnum = it.colnum;
					return *this;
			}

			// Comparison
			BOOST_UBLAS_INLINE
				bool operator == (const const_iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return((rownum==it.rownum)&&(colnum==it.colnum));
			}
			BOOST_UBLAS_INLINE
				bool operator < (const const_iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return(rownum < it.rownum);
			}

		private:
			size_type rownum;
			size_type colnum;
		};
#endif

		BOOST_UBLAS_INLINE
			const_iterator1 begin1 () const {
				return find1 (0, 0, 0);
		}
		BOOST_UBLAS_INLINE
			const_iterator1 end1 () const {
				return find1 (0, size1 (), 0);
		}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
		class iterator1:
			public container_reference<vector_tile>,
			public random_access_iterator_base<typename iterator_restrict_traits<
			typename subiterator1_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
			iterator1, value_type> {
		public:
			typedef typename subiterator1_type::value_type value_type;
			typedef typename subiterator1_type::difference_type difference_type;
			typedef typename subiterator1_type::reference reference;
			typedef typename subiterator1_type::pointer pointer;

			typedef iterator2 dual_iterator_type;
			typedef reverse_iterator2 dual_reverse_iterator_type;

			// Construction and destruction
			BOOST_UBLAS_INLINE
				iterator1 ():
			container_reference<self_type> (), rownum (), colnum () {}
			BOOST_UBLAS_INLINE
				iterator1 (self_type &m, size_type row_number, size_type column_number):
			container_reference<self_type> (m), rownum(row_number), colnum(column_number) {}

			// Arithmetic
			BOOST_UBLAS_INLINE
				iterator1 &operator ++ () {
					++ rownum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator1 &operator -- () {
					-- rownum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator1 &operator += (difference_type n) {
					rownum += n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator1 &operator -= (difference_type n) {
					rownum -= n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				difference_type operator - (const iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (rownum - it.rownum);
			}

			// Dereference
			BOOST_UBLAS_INLINE
				reference operator * () const {
					size_type i = index1 ();
					size_type j = index2 ();
					BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
					BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
					return (*this) () (i, j);
			}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				iterator2 begin () const {
					return (*this) ().find2 (1, index1 (), 0);
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				iterator2 end () const {
					return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				reverse_iterator2 rbegin () const {
					return reverse_iterator2 (end ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				reverse_iterator2 rend () const {
					return reverse_iterator2 (begin ());
			}
#endif

			// Indices
			BOOST_UBLAS_INLINE
				size_type index1 () const {
					return rownum;
			}
			BOOST_UBLAS_INLINE
				size_type index2 () const {
					return colnum;
			}

			// Assignment
			BOOST_UBLAS_INLINE
				iterator1 &operator = (const iterator1 &it) {
					container_reference<self_type>::assign (&it ());
					rownum = it.rownum;
					colnum = it.colnum;
					return *this;
			}

			// Comparison
			BOOST_UBLAS_INLINE
				bool operator == (const iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return ((rownum==it.rownum)&&(colnum==it.colnum));
			}
			BOOST_UBLAS_INLINE
				bool operator < (const iterator1 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (rownum - it.rownum);
			}

		private:
			size_type rownum;
			size_type colnum; 
			friend class const_iterator1;
		};
#endif

		BOOST_UBLAS_INLINE
			iterator1 begin1 () {
				return find1 (0, 0, 0);
		}
		BOOST_UBLAS_INLINE
			iterator1 end1 () {
				return find1 (0, size1 (), 0);
		}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
		class const_iterator2:
			public container_const_reference<vector_tile>,
			public random_access_iterator_base<packed_random_access_iterator_tag,
			const_iterator2, value_type> {
		public:
			typedef typename iterator_restrict_traits<typename const_subiterator2_type::iterator_category,
				packed_random_access_iterator_tag>::iterator_category iterator_category;
			typedef typename const_subiterator2_type::value_type value_type;
			typedef typename const_subiterator2_type::difference_type difference_type;
			typedef typename const_subiterator2_type::reference reference;
			typedef typename const_subiterator2_type::pointer pointer;

			typedef const_iterator1 dual_iterator_type;
			typedef const_reverse_iterator1 dual_reverse_iterator_type;

			// Construction and destruction
			BOOST_UBLAS_INLINE
				const_iterator2 ():
			container_const_reference<self_type> (), rownum(), colnum () {}
			BOOST_UBLAS_INLINE
				const_iterator2 (const self_type &m, const size_type row_number, const size_type column_number):
			container_const_reference<self_type> (m), rownum (row_number), colnum (column_number) {}
			BOOST_UBLAS_INLINE
				const_iterator2 (const iterator2 &it):
			container_const_reference<self_type> (it ()), rownum (it.rownum), colnum (it.colnum) {}

			// Arithmetic
			BOOST_UBLAS_INLINE
				const_iterator2 &operator ++ () {
					++ colnum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator2 &operator -- () {
					-- colnum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator2 &operator += (difference_type n) {
					colnum += n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				const_iterator2 &operator -= (difference_type n) {
					colnum -= n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				difference_type operator - (const const_iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (colnum - it.colnum);
			}

			// Dereference
			BOOST_UBLAS_INLINE
				const_reference operator * () const {
					size_type i = index1 ();
					size_type j = index2 ();
					BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
					BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
					return (*this) () (i, j);
			}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_iterator1 begin () const {
					return (*this) ().find1 (1, 0, index2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_iterator1 end () const {
					return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_reverse_iterator1 rbegin () const {
					return const_reverse_iterator1 (end ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				const_reverse_iterator1 rend () const {
					return const_reverse_iterator1 (begin ());
			}
#endif

			// Indices
			BOOST_UBLAS_INLINE
				size_type index1 () const {
					return rownum;
			}
			BOOST_UBLAS_INLINE
				size_type index2 () const {
					return colnum;
			}

			// Assignment
			BOOST_UBLAS_INLINE
				const_iterator2 &operator = (const const_iterator2 &it) {
					container_const_reference<self_type>::assign (&it ());
					rownum = it.rownum;
					colnum = it.colnum;
					return *this;
			}

			// Comparison
			BOOST_UBLAS_INLINE
				bool operator == (const const_iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return ((rownum==it.rownum)&&(colnum==it.colnum));
			}
			BOOST_UBLAS_INLINE
				bool operator < (const const_iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (colnum < it.colnum);
			}

		private:
			size_type rownum;
			size_type colnum;
		};
#endif

		BOOST_UBLAS_INLINE
			const_iterator2 begin2 () const {
				return find2 (0, 0, 0);
		}
		BOOST_UBLAS_INLINE
			const_iterator2 end2 () const {
				return find2 (0, 0, size2 ());
		}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
		class iterator2:
			public container_reference<vector_tile>,
			public random_access_iterator_base<typename iterator_restrict_traits<
			typename subiterator2_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
			iterator2, value_type> {
		public:
			typedef typename subiterator2_type::value_type value_type;
			typedef typename subiterator2_type::difference_type difference_type;
			typedef typename subiterator2_type::reference reference;
			typedef typename subiterator2_type::pointer pointer;

			typedef iterator1 dual_iterator_type;
			typedef reverse_iterator1 dual_reverse_iterator_type;

			// Construction and destruction
			BOOST_UBLAS_INLINE
				iterator2 ():
			container_reference<self_type> (), rownum () , colnum(){}
			BOOST_UBLAS_INLINE
				iterator2 (self_type &m, size_type row_number, size_type column_number):
			container_reference<self_type> (m), rownum(row_number), colnum (column_number)  {}

			// Arithmetic
			BOOST_UBLAS_INLINE
				iterator2 &operator ++ () {
					++ colnum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator2 &operator -- () {
					-- colnum;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator2 &operator += (difference_type n) {
					colnum += n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				iterator2 &operator -= (difference_type n) {
					colnum -= n;
					return *this;
			}
			BOOST_UBLAS_INLINE
				difference_type operator - (const iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return colnum - it.colnum;
			}

			// Dereference
			BOOST_UBLAS_INLINE
				reference operator * () const {
					size_type i = index1 ();
					size_type j = index2 ();
					return (*this) () (i,j);
			}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				iterator1 begin () const {
					return (*this) ().find1 (1, 0, index2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				iterator1 end () const {
					return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				reverse_iterator1 rbegin () const {
					return reverse_iterator1 (end ());
			}
			BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
				typename self_type::
#endif
				reverse_iterator1 rend () const {
					return reverse_iterator1 (begin ());
			}
#endif

			// Indices
			BOOST_UBLAS_INLINE
				size_type index1 () const {
					return rownum;
			}
			BOOST_UBLAS_INLINE
				size_type index2 () const {
					return colnum;
			}

			// Assignment
			BOOST_UBLAS_INLINE
				iterator2 &operator = (const iterator2 &it) {
					container_reference<self_type>::assign (&it ());
					rownum = it.rownum;
					colnum = it.colnum;
					return *this;
			}

			// Comparison
			BOOST_UBLAS_INLINE
				bool operator == (const iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return ((rownum==it.rownum)&&(colnum==it.colnum));
			}
			BOOST_UBLAS_INLINE
				bool operator < (const iterator2 &it) const {
					BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
					return (colnum < it.colnum);
			}

		private:
			size_type rownum;
			size_type colnum;
			friend class const_iterator2;
		};
#endif

		BOOST_UBLAS_INLINE
			iterator2 begin2 () {
				return find2 (0, 0, 0);
		}
		BOOST_UBLAS_INLINE
			iterator2 end2 () {
				return find2 (0, 0, size2 ());
		}

		// Reverse iterators

		BOOST_UBLAS_INLINE
			const_reverse_iterator1 rbegin1 () const {
				return const_reverse_iterator1 (end1 ());
		}
		BOOST_UBLAS_INLINE
			const_reverse_iterator1 rend1 () const {
				return const_reverse_iterator1 (begin1 ());
		}

		BOOST_UBLAS_INLINE
			reverse_iterator1 rbegin1 () {
				return reverse_iterator1 (end1 ());
		}
		BOOST_UBLAS_INLINE
			reverse_iterator1 rend1 () {
				return reverse_iterator1 (begin1 ());
		}

		BOOST_UBLAS_INLINE
			const_reverse_iterator2 rbegin2 () const {
				return const_reverse_iterator2 (end2 ());
		}
		BOOST_UBLAS_INLINE
			const_reverse_iterator2 rend2 () const {
				return const_reverse_iterator2 (begin2 ());
		}

		BOOST_UBLAS_INLINE
			reverse_iterator2 rbegin2 () {
				return reverse_iterator2 (end2 ());
		}
		BOOST_UBLAS_INLINE
			reverse_iterator2 rend2 () {
				return reverse_iterator2 (begin2 ());
		}

	private:
		matrix_closure_type data_;
		size_type num_rows;
		size_type num_cols;
	};


	template<class V>
	BOOST_UBLAS_INLINE
		vector_tile<V> matvec(V& data) {
			return vector_tile<V>(data);
	};

	template<class V>
	BOOST_UBLAS_INLINE
		vector_tile<V> matvec(V& data,int r,int c){
			return vector_tile<V>(data,r,c);
	};

	template<class V>
	BOOST_UBLAS_INLINE
		const vector_tile<const V> matvec(const V& data) {
			return vector_tile<const V>(data);
	};


}}}

#endif //_BOOST_TILE_MATRIX_FROM_VECTOR_
