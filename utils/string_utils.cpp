#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>

namespace mine{
	std::string int2string(unsigned int n);
}

	std::string mine::int2string(unsigned int n)
	{
		std::string s("");
		int remainder;
		do{
			remainder = n % 10;
			n /= 10;
			s += (char)(remainder + 48);
		}while(n > 0);
		std::string reversed(s.rbegin(), s.rend());
		return reversed;
	}

#endif