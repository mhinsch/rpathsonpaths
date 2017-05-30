#ifndef UTIL_H
#define UTIL_H

#include <stdexcept>

template <class STR>
inline void assert_failed(const STR & str)
	{
	throw std::runtime_error(str);
	}

#ifndef NO_ASSERT

#define ASSERT_STRFY(arg) #arg
#define ASSERT_TOSTR(arg) ASSERT_STRFY(arg)

#define myassert(cond) ((cond) ? (void)0 : assert_failed(__FILE__ ":" ASSERT_TOSTR(__LINE__) " Condition '"  #cond "' failed!"))

#else

#define myassert(cond) ((void)0)

#endif

#endif	// UTIL_H
