#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"


using namespace tbb;
using namespace std;

#define THRESHOLD 128
#define FILTER_SIZE				5
#define DIM_REDUCTION			9
#define CUTOFF					350	

class ParallelPrewitt {
	int* input;
	int* output;
	int height;
	int width;
	int* filterHor;
	int* filterVer;
public:
	ParallelPrewitt(int* input_, int* output_, int width_, int height_, int* filterHor_, int* filterVer_) :
		input(input_), output(output_), width(width_), height(height_) 
	{
		filterHor = filterHor_;
		filterVer = filterVer_;
	}

	void operator() (const blocked_range<int>& r) const {
		int* a = output;
		for (int i = r.begin(); i != r.end(); ++i) {
			for (int j = DIM_REDUCTION; j < height - DIM_REDUCTION; j++) {
				int Fx = 0;
				int Fy = 0;
				int F = 0;
				for (int m = 0; m < FILTER_SIZE; m++) {
					for (int n = 0; n < FILTER_SIZE; n++) {
						Fx += input[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterHor[m * FILTER_SIZE + n];
						Fy += input[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterVer[m * FILTER_SIZE + n];
					}
				}
				F = abs(Fx) + abs(Fy);

				if (F > THRESHOLD) {
					output[j * width + i] = 255;
				}
				else {
					output[j * width + i] = 0;
				}
			}
		}
	}
};