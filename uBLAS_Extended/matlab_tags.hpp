#ifndef UBLAS_MATLAB_TAGS_H
#define UBLAS_MATLAB_TAGS_H

namespace boost { namespace numeric { namespace ublas { 
	template <class D = std::ptrdiff_t>
	struct last_elem {
		typedef D  difference_type;

		BOOST_UBLAS_INLINE
			last_elem():offset_(0){}

		explicit 	BOOST_UBLAS_INLINE
			last_elem(difference_type offset):offset_(offset){}

		BOOST_UBLAS_INLINE
			last_elem(const last_elem<D>& other):offset_(other.offset_){}

		BOOST_UBLAS_INLINE
			difference_type get_offset()const{
				return offset_;
		}
	private:
		difference_type offset_;
	};

	template <class T>
	BOOST_UBLAS_INLINE
		last_elem<T> operator + (const last_elem<T>& elem, T offset){
			return last_elem<T>(offset + elem.get_offset());
	}

	template <class T>
	BOOST_UBLAS_INLINE
		last_elem<T> operator + (T offset, const last_elem<T>& elem){
			return last_elem<T>(offset + elem.get_offset());
	}

	template <class T>
	BOOST_UBLAS_INLINE
		last_elem<T> operator - (const last_elem<T>& elem, T offset){
			return last_elem<T>(elem.get_offset() - offset);
	}

	template <class T>
	BOOST_UBLAS_INLINE
		last_elem<T> operator - (T offset, const last_elem<T>& elem){
			return last_elem<T>(offset - elem.get_offset());
	}
	struct all_elem{};
	struct rev_all{};

    template <class Z = std::size_t, class D = std::ptrdiff_t>
	class basic_range_ext {
	public:
		typedef Z size_type;
		typedef D difference_type;
		basic_range_ext(const last_elem<D>& start, const last_elem<D>& end):
							start_(start.get_offset()),
							end_(end.get_offset()),
							bIsStartOffset_(true),
							bIsEndOffset_(true){}
		basic_range_ext(const last_elem<D>& start, size_type end):
							start_(start.get_offset()),
							end_((difference_type)end),
							bIsStartOffset_(true),
							bIsEndOffset_(false){}	
		basic_range_ext(size_type start, const last_elem<D>& end):
							start_((difference_type)start),
							end_(end.get_offset()),
							bIsStartOffset_(false),
							bIsEndOffset_(true){}
	private:
		difference_type start_;
		difference_type end_;
		bool bIsStartOffset_;
		bool bIsEndOffset_;

	};

    template <class Z = std::size_t, class D = std::ptrdiff_t>
	class basic_slice_ext {
	public:
		typedef Z size_type;
		typedef D difference_type;
		basic_slice_ext(const last_elem<D>& start, difference_type stride, const last_elem<D>& end):
							start_(start.get_offset()),
							end_(end.get_offset()),
							bIsStartOffset_(true),
							bIsEndOffset_(true),
							stride_(stride){}
		basic_slice_ext(const last_elem<D>& start, difference_type stride, size_type end):
							start_(start.get_offset()),
							end_((difference_type)end),
							bIsStartOffset_(true),
							bIsEndOffset_(false),
							stride_(stride){}	
		basic_slice_ext(size_type start, difference_type stride, const last_elem<D>& end):
							start_((difference_type)start),
							end_(end.get_offset()),
							bIsStartOffset_(false),
							bIsEndOffset_(true),
							stride_(stride){}
		basic_slice_ext(const last_elem<D>& start, const last_elem<D>& end):
							start_(start.get_offset()),
							end_(end.get_offset()),
							bIsStartOffset_(true),
							bIsEndOffset_(true),
							stride_(1){}
		basic_slice_ext(size_type start, difference_type stride, size_type end):
							start_((difference_type)start),
							end_((difference_type)end),
							bIsStartOffset_(false),
							bIsEndOffset_(false),
							stride_(stride){}	
		basic_slice_ext(const last_elem<D>& start, size_type end):
							start_(start.get_offset()),
							end_((difference_type)end),
							bIsStartOffset_(true),
							bIsEndOffset_(false),
							stride_(1){}	
		basic_slice_ext(size_type start, const last_elem<D>& end):
							start_((difference_type)start),
							end_(end.get_offset()),
							bIsStartOffset_(false),
							bIsEndOffset_(true),
							stride_(1){}
		basic_slice_ext(size_type start, size_type end):
							start_((difference_type)start),
							end_((difference_type)end),
							bIsStartOffset_(false),
							bIsEndOffset_(false),
							stride_(1){}	
	public:
		difference_type start_;
		difference_type end_;
		difference_type stride_;
		bool bIsStartOffset_;
		bool bIsEndOffset_;

	};
}}}
#endif //UBLAS_MATLAB_TAGS_H
