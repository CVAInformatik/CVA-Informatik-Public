#pragma once

#ifdef OS_WINDOWS    // windows
#define WIN
typedef unsigned __int64  u64;
typedef __int64  s64;
#else   //linux
#define NOTWIN
#include <stdlib.h>
#include <cstdint>
typedef  uint64_t  u64;
typedef  int64_t  s64;
#endif


class GrayCodeGeneratorClass
{
public:
	GrayCodeGeneratorClass() {};
	~GrayCodeGeneratorClass() {};

	inline void Init(unsigned int _width, u64 FirstValue = 0) {
		
		if (0 < _width && _width < 65) {
			switch (_width)
			{
			case 64: mask = 0xFFFFFFFFFFFFFFFF;
				break;
			default: mask = ((u64)1) << _width;
				mask--;
				break;
			}

			u64 temp = (FirstValue & mask) >> 1;
			state = (FirstValue & mask);
			while (temp) { 	state ^= temp; temp = temp >> 1; }

		}
		else { mask = 0; state = 0; }
	};

	inline u64 Next() { u64 res = state ^ (state >> 1) ; state = (state + 1) & mask; return mask & res;	};

private:
	u64 state;
	u64 mask;
};
