#ifndef UBLAS_FIND_EXT_H
#define UBLAS_FIND_EXT_H

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace boost { namespace numeric { namespace ublas {

	// find(boolean vector expression exprCond, vector expression for return values exprSearch)  
	// returns exprSearch(elements where exprCond == true)
	template<class E, class E1>
	vector_indirect<const E1, indirect_array<unbounded_array<typename E1::size_type> > >
		find(const vector_expression<E>& exprCond, 
			const vector_expression<E1>& exprSearch, 
			typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef typename E1::size_type size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type;
			vector<size_type> index;
			return exprSearch()(find(exprCond, index));
	}
	// find(boolean vector expression exprCond, 
    //				vector expression for return values exprSearch, indices where boolean expression evaluates to true)
	// returns exprSearch(elements where exprCond == true)
	template<class E, class E1>
	vector_indirect<const E1, indirect_array<unbounded_array<typename E1::size_type> > >
		find(const vector_expression<E>& exprCond, 
			const vector_expression<E1>& exprSearch, 
			vector<typename E1::size_type>& index, 
			typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef typename E1::size_type size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type;
			return exprSearch()(find(exprCond, index));
	}

	// find(boolean vector expression exprCond, indices where boolean expression evaluates to true)
	// returns an indirect array that can be used to index 
	template<class E, class S>
	indirect_array<unbounded_array<typename E::size_type> > 
		find(const vector_expression<E>& exprCond, 
			vector<S>& index, 
			typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef S size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type;
			typedef vector_indirect<const E, indirect_type> vector_indirect_type;
			size_type count = 0;
			size_type nMaxFind = exprCond().size();
			index.resize(nMaxFind, false);
			for(size_type ii = 0; ii < nMaxFind; ii++){
				if(exprCond()(ii)) 
					index(count++) = ii;
 			}
			indirect_type indIndex(count, index.data());
			return indIndex;
	}

	// find(boolean matrix expression exprCond, matrix expression for return values exprSearch)  
	// returns exprSearch(elements where exprCond == true)
	template<class E, class E1>
	matrix_vector_indirect<const E1, indirect_array<unbounded_array<typename E1::size_type> > >
		find(const matrix_expression<E>& exprCond, 
			 const matrix_expression<E1>& exprSearch, 
			typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef typename E1::size_type size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type; 
			typedef matrix_vector_indirect<const E1,  indirect_type>  matrix_indirect_type;
			vector<size_type> index1, index2;
			find(exprCond, index1, index2);
			return matrix_indirect_type(exprSearch(), 
				indirect_type(index1.size(), index1.data()), indirect_type(index2.size(), index2.data()));
	}

	// find(boolean matrix expression exprCond, 
    //				matrix expression for return values exprSearch, 
	//				row indices where boolean expression evaluates to true,
	//				column indices where boolean expression evaluates to true,)
	// returns exprSearch(elements where exprCond == true) - linear indexed matrix
	template<class E, class E1>
	matrix_vector_indirect<const E1, indirect_array<unbounded_array<typename E1::size_type> > >
		find(const matrix_expression<E>& exprCond, 
			 const matrix_expression<E1>& exprSearch, 
			 vector<typename E1::size_type>& index1,
			 vector<typename E1::size_type>& index2,
			typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef typename E1::size_type size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type; 
			typedef matrix_vector_indirect<const E1,  indirect_type>  matrix_indirect_type;
			find(exprCond, index1, index2);
			return matrix_indirect_type(exprSearch(), 
				indirect_type(index1.size(), index1.data()), indirect_type(index2.size(), index2.data()));
	}

	// find(boolean matrix expression exprCond, 
	//				row indices where boolean expression evaluates to true,
	//				column indices where boolean expression evaluates to true)
	template<class E, class S>
	void find(const matrix_expression<E>& expr, 
			 vector<S>& index1,
			 vector<S>& index2,
			 typename boost::enable_if<typename boost::is_same<typename E::value_type, bool>::type>::type *dummy = 0){
			typedef S size_type;
			typedef indirect_array<unbounded_array<size_type> > indirect_type;
			typedef matrix_vector_indirect<const E, indirect_type> vector_indirect_type;
			size_type count = 0, nSize1 = expr().size1(), nSize2 = expr().size2();
			size_type nMaxFind = nSize1*nSize2;
			index1.resize(nMaxFind, false);
			index2.resize(nMaxFind, false);
			for(size_type ii = 0; ii < nSize1; ii++)
				for(size_type jj = 0; jj < nSize2; jj++){
					if(expr()(ii, jj)) {
						index1(count) = ii;
						index2(count++) = jj;
					}
 			}
			index1.resize(count);
			index2.resize(count);
	}

}}};
#endif //UBLAS_FIND_EXT_H
