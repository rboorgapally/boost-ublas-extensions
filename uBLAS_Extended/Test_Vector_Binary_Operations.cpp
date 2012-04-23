	#include "stdafx.h"
	#include "vector_expression_ext.hpp"
	#include <boost/numeric/ublas/vector.hpp>
	#include <boost/numeric/ublas/io.hpp>
	void testVectorBinaryOperations () {
	  using namespace boost::numeric::ublas;
	  vector<double> v1(10), v2(10);
	  vector<bool> b1(10), b2(10);
	  for (unsigned i = 0; i < v1.size (); ++ i) {
		v1(i) = 1;
		v2(i) = 1;
		b1 (i) = true;
		b2 (i) = false;
	  }
	  std::cout << (2 + v1) << std::endl;
	  std::cout << (v1 + 4) << std::endl;
	  std::cout << (2 - v1) << std::endl;
	  std::cout << v1 - 4 << std:: endl;
	  std::cout << (v1 ^ 2) << std:: endl;
	  std::cout << (v1 > 1) << std:: endl;
	  std::cout << (1 > v1) << std::endl;
	  std::cout << (v1 > v2) << std::endl;
	  std::cout << (v1 < 1) << std::endl;
	  std::cout << (1 < v1) << std::endl;
	  std::cout << (v1 < v2) << std::endl;
	  std::cout << (v1 >= 1) << std::endl;
	  std::cout << (1 >= v1) << std::endl;
	  std::cout << (v1 >= v2) << std::endl;
	  std::cout << (v1 <= 1) << std::endl;
	  std::cout << (1 <= v1) << std::endl;
	  std::cout << (v1 <= v1) << std::endl;
	  std::cout << (v1 != 1) << std::endl;
	  std::cout << (1 != v1) << std::endl;
	  std::cout << (v2 != v1) << std::endl;
	  std::cout << (b1 && true) << std::endl;
	  std::cout << (true && b1) << std::endl;
	  std::cout << (b1 && b2) << std::endl;
	  std::cout << (b2 || true) << std::endl;
	  std::cout << (true || b2) << std::endl;
	  std::cout << (b1 || b2) << std::endl;
	  std::cout << (!b2) << std::endl;

	  std::cout << "Test for Vector binary operators complete" << std::endl;
	}
