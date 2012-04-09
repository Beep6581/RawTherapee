#ifndef _MYMATH_
#define _MYMATH_
#include <cmath>
#include <algorithm>

#include "rtengine.h"

namespace rtengine {
	template <typename _Tp>
	inline const _Tp SQR (_Tp x) {
//		return std::pow(x,2); Slower than:
		return (x*x);
	}

	template<typename _Tp>
	inline const _Tp& min(const _Tp& a, const _Tp& b) {
		return std::min(a,b);
	}

	template<typename _Tp>
	inline const _Tp& max(const _Tp& a, const _Tp& b) {
		return std::max(a,b);
	}

	template<typename _Tp>
	inline const _Tp& min(const _Tp& a, const _Tp& b, const _Tp& c) {
		return std::min(c,std::min(a,b));
	}

	template<typename _Tp>
	inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c) {
		return std::max(c,std::max(a,b));
	}

	template<typename _Tp>
	inline const _Tp& min(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d) {
		return std::min(d,std::min(c,std::min(a,b)));
	}

	template<typename _Tp>
	inline const _Tp& max(const _Tp& a, const _Tp& b, const _Tp& c, const _Tp& d) {
		return std::max(d,std::max(c,std::max(a,b)));
	}
}
#endif
