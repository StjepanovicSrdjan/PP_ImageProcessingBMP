#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "BitmapRawConverter.h"
#include "tbb/tick_count.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include <tbb/task_group.h>
#include "ParallelPrewitt.cpp"

#define __ARG_NUM__				6
#define FILTER_SIZE				3
#define DIM_REDUCTION			1
#define THRESHOLD				128
#define CUTOFF					250	
#define DEPTH					1

using namespace std;
using namespace tbb;

ofstream output_prewitt_fs;
ofstream output_prewitt_co;
ofstream output_edge_co;

// Prewitt operators
int filterHor[FILTER_SIZE * FILTER_SIZE] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
int filterVer[FILTER_SIZE * FILTER_SIZE] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};

//int filterHor[FILTER_SIZE * FILTER_SIZE] = { 9, 9, 9, 9, 9, 9, 5, 5, 5, 9, -7, -3, 0, -3, -7, -7, -3, -3, -3, -7, -7, -7, -7, -7, -7 };
//int filterVer[FILTER_SIZE * FILTER_SIZE] = { 9, 9, -7, -7, -7, 9, 5, -3, -3, -7, 9, 5, 0, -3, -7, 9, 5, -3, -3, -7, 9, 9, -7, -7, -7 };

//int filterHor[FILTER_SIZE * FILTER_SIZE] = { 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1 };
//int filterVer[FILTER_SIZE * FILTER_SIZE] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };

//int filterHor[FILTER_SIZE * FILTER_SIZE] = { 1, 1, 1, 0, -1, -1, -1, 1, 1, 1, 0, -1, -1, -1, 1, 1, 1, 0, -1, -1, -1, 1, 1, 1, 0, -1, -1, -1 , 1, 1, 1, 0, -1, -1, -1 , 1, 1, 1, 0, -1, -1, -1 , 1, 1, 1, 0, -1, -1, -1 };
//int filterVer[FILTER_SIZE * FILTER_SIZE] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

/**
* @brief Serial version of edge detection algorithm implementation using Prewitt operator
* @param inBuffer buffer of input image
* @param outBuffer buffer of output image
* @param width image width
* @param height image height
*/
void filter_serial_prewitt(int *inBuffer, int *outBuffer, int width, int height)  //TODO obrisati
{
	for (int i = DIM_REDUCTION; i < width - DIM_REDUCTION; i++) {
		for (int j = DIM_REDUCTION; j < height - DIM_REDUCTION; j++) {
			int Fx = 0;
			int Fy = 0;
			int F = 0;

			for (int m = 0; m < FILTER_SIZE; m++) {
				for (int n = 0; n < FILTER_SIZE; n++) {
					Fx += inBuffer[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterHor[m * FILTER_SIZE + n];
					Fy += inBuffer[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterVer[m * FILTER_SIZE + n];
				}
			}
			F = abs(Fx) + abs(Fy);
			
			if (F > THRESHOLD) {
				outBuffer[j * width + i] = 255;
			}
			else {
				outBuffer[j * width + i] = 0;
			}
		}
	}


}
/**
* @brief parallel version of prewitt algorithm using parallel_for 
* @param inBuffer buffer of input image* @param outBuffer buffer of output image
* @param width image width
* @param height image height
*/
void filter_prewitt_parallelFor(int* inBuffer, int* outBuffer, int width, int height) {
	ParallelPrewitt pp(inBuffer, outBuffer, width, height, (int*) filterHor, (int*) filterVer);
	parallel_for(blocked_range<int>(DIM_REDUCTION, width - DIM_REDUCTION), pp);
}


/**
* @brief Parallel version of edge detection algorithm implementation using Prewitt operator
* 
* @param inBuffer buffer of input image
* @param outBuffer buffer of output image
* @param width image width
* @param height image height
*/
void filter_parallel_prewitt(int *inBuffer, int *outBuffer, int width, int height, int xStart, int xEnd, int yStart, int yEnd)
{
	

	task_group g;
	if (xEnd - xStart <= CUTOFF) {
		for (int i = xStart; i < xEnd; i++) {
			for (int j = yStart; j < yEnd; j++) {
				int Fx = 0;
				int Fy = 0;
				int F = 0;
				
				for (int m = 0; m < FILTER_SIZE; m++) {
					for (int n = 0; n < FILTER_SIZE; n++) {
						Fx += inBuffer[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterHor[m * FILTER_SIZE + n];
						Fy += inBuffer[(j + filterHor[m * FILTER_SIZE + n]) * width + (i + filterVer[m * FILTER_SIZE + n])] * filterVer[m * FILTER_SIZE + n];
					}
				}
				F = abs(Fx) + abs(Fy);

				if (F > THRESHOLD) {
					outBuffer[j * width + i] = 255;
				}
				else {
					outBuffer[j * width + i] = 0;
				}
			}
		}
	}
	else {
		g.run([&] {filter_parallel_prewitt(inBuffer, outBuffer, width, height, xStart, (xStart + xEnd) / 2, yStart, (yStart + yEnd) / 2); });
		g.run([&] {filter_parallel_prewitt(inBuffer, outBuffer, width, height, xStart, (xStart + xEnd) / 2, (yStart + yEnd) / 2, yEnd); });
		g.run([&] {filter_parallel_prewitt(inBuffer, outBuffer, width, height, (xStart + xEnd) / 2, xEnd, yStart, (yStart + yEnd) / 2); });
		g.run([&] {filter_parallel_prewitt(inBuffer, outBuffer, width, height, (xStart + xEnd) / 2, xEnd, (yStart + yEnd) / 2, yEnd); });
		g.wait();
	}


}


void  thresholdCutOff(int* inBuffer, int width, int height)
{
	for (int i = DEPTH; i < width - DEPTH; i++) {
		for (int j = 1; j < height - 1; j++) {
			inBuffer[j * width + i] = (inBuffer[j * width + i] < THRESHOLD) ? 0 : 1;
		}
	}
}


/**
* @brief Serial version of edge detection algorithm
* @param inBuffer buffer of input image
* @param outBuffer buffer of output image
* @param width image width
* @param height image height
*/
void filter_serial_edge_detection(int *inBuffer, int *outBuffer, int width, int height)	//TODO obrisati
{
	const int up = -width;
	const int down = width;
	const int right = 1;
	const int left = -1;
	const int upRight = 1 - width;
	const int upLeft = -1 - width;
	const int downRight = 1 + width;
	const int downLeft = width - 1;

	int currentIndex;
	for (int i = DEPTH; i < width - DEPTH; i++) {
		for (int j = DEPTH; j < height - DEPTH; j++) {
			int p = 0;
			int o = 0;
			currentIndex = j * width + i;
			for (int n = 1; n <= DEPTH; n++) {
				if (inBuffer[currentIndex + n*up] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*up] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*down] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*down] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*right] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*right] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*left] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*left] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*upRight] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*upRight] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*upLeft] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*upLeft] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*downLeft] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*downLeft] == 0) {
					o++;
				}

				if (inBuffer[currentIndex + n*downRight] == 1) {
					p++;
				}
				else if (inBuffer[currentIndex + n*downRight] == 0) {
					o++;
				}
			}
			if (p != 0) p = 1;
			if (o != 0) o = 1;

			outBuffer[currentIndex] = (abs(p) - abs(o) == 0) ? 255 : 0;
		}

	}

}

/**
* @brief Parallel version of edge detection algorithm
* 
* @param inBuffer buffer of input image
* @param outBuffer buffer of output image
* @param width image width
* @param height image height
* @param xStart x axis start index
* @param xEnd x axis end index
* @param yStart y axis start index
* @param yEnd y axis end index
*/
void filter_parallel_edge_detection(int* inBuffer, int* outBuffer, int width, int height, int xStart, int xEnd, int yStart, int yEnd)
{
	task_group g;
	if (xEnd - xStart <= CUTOFF) {
		const int up = -width;
		const int down = width;
		const int right = 1;
		const int left = -1;
		const int upRight = 1 - width;
		const int upLeft = -1 - width;
		const int downRight = 1 + width;
		const int downLeft = width - 1;

		int currentIndex;
		for (int i = xStart; i < xEnd; i++) {
			for (int j = yStart; j < yEnd; j++) {
				int p = 0;
				int o = 0;
				currentIndex = j * width + i;
				for (int n = 1; n <= DEPTH; n++) {
					if (inBuffer[currentIndex + n * up] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * up] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * down] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * down] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * right] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * right] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * left] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * left] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * upRight] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * upRight] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * upLeft] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * upLeft] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * downLeft] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * downLeft] == 0) {
						o++;
					}

					if (inBuffer[currentIndex + n * downRight] == 1) {
						p++;
					}
					else if (inBuffer[currentIndex + n * downRight] == 0) {
						o++;
					}
				}
				if (p != 0) p = 1;
				if (o != 0) o = 1;

				outBuffer[currentIndex] = (abs(p) - abs(o) == 0) ? 255 : 0;
			}
		}
	}
	else {
		g.run([&] {filter_parallel_edge_detection(inBuffer, outBuffer, width, height, xStart, (xStart + xEnd) / 2, yStart, (yStart + yEnd) / 2); });
		g.run([&] {filter_parallel_edge_detection(inBuffer, outBuffer, width, height, xStart, (xStart + xEnd) / 2, (yStart + yEnd) / 2, yEnd); });
		g.run([&] {filter_parallel_edge_detection(inBuffer, outBuffer, width, height, (xStart + xEnd) / 2, xEnd, yStart, (yStart + yEnd) / 2); });
		g.run([&] {filter_parallel_edge_detection(inBuffer, outBuffer, width, height, (xStart + xEnd) / 2, xEnd, (yStart + yEnd) / 2, yEnd); });
		g.wait();
	}
}



/**
* @brief Function for running test.
*
* @param testNr test identification, 1: for serial version, 2: for parallel version
* @param ioFile input/output file, firstly it's holding buffer from input image and than to hold filtered data
* @param outFileName output file name
* @param outBuffer buffer of output image
* @param width image width
* @param height image height
*/


void run_test_nr(int testNr, BitmapRawConverter* ioFile, char* outFileName, int* outBuffer, unsigned int width, unsigned int height, ofstream& file1, ofstream& file2)
{

	// TODO: start measure
	string output = "";
	tick_count start = tick_count::now();

	switch (testNr)
	{
		case 1:
			cout << "Running serial version of edge detection using Prewitt operator" << endl;
			filter_serial_prewitt(ioFile->getBuffer(), outBuffer, width, height);
			output += "Serial : ";
			break;
		case 2:
			cout << "Running parallel version of edge detection using Prewitt operator" << endl;
			filter_parallel_prewitt(ioFile->getBuffer(), outBuffer, width, height, DIM_REDUCTION, width - DIM_REDUCTION, DIM_REDUCTION, height - DIM_REDUCTION);
			output += "Parallel : ";
			break;
		case 3:
			cout << "Running serial version of edge detection" << endl;
			filter_serial_edge_detection(ioFile->getBuffer(), outBuffer, width, height);
			output += "Serial : ";
			break;
		case 4:
			cout << "Running parallel version of edge detection" << endl;
			filter_parallel_edge_detection(ioFile->getBuffer(), outBuffer, width, height, DEPTH, width - DEPTH, DEPTH, height - DEPTH);
			output += "Parallel : ";
			break;
		case 5:
			cout << "Running parallel version of edge detection using Prewitt operator, parallelized using parallel_for" << endl;
			filter_prewitt_parallelFor(ioFile->getBuffer(), outBuffer, width, height);
			output += "ParallelFor : ";
			break;
		default:
			cout << "ERROR: invalid test case, must be 1, 2, 3 or 4!";
			break;
	}
	// TODO: end measure and display time

	tick_count end = tick_count::now();
	output += to_string((end - start).seconds());
	cout << "\n" << testNr << " : " << (end - start).seconds() << endl;

	if (testNr != 3 && testNr != 4) {
		file1 << output << "\n";
		file2 << output << "\n";
	}
	else {
		file1 << output << "\n";
	}

	ioFile->setBuffer(outBuffer);
	ioFile->pixelsToBitmap(outFileName);
}

/**
* @brief Print program usage.
*/
void usage()
{
	cout << "\n\ERROR: call program like: " << endl << endl; 
	cout << "ProjekatPP.exe";
	cout << " input.bmp";
	cout << " outputSerialPrewitt.bmp";
	cout << " outputParallelPrewitt.bmp";
	cout << " outputSerialEdge.bmp";
	cout << " outputParallelEdge.bmp" << endl << endl;
}


int main(int argc, char * argv[])
{
	string file_name = "cutOff_" + to_string(CUTOFF) + "-filterSize_" + to_string(FILTER_SIZE) + ".txt";
	cout << file_name;
	output_prewitt_fs.open("../test_data/PrewittDiffFilterSize/cutOff_" + to_string(CUTOFF) + "-filterSize_" + to_string(FILTER_SIZE) + ".txt");
	output_prewitt_co.open("../test_data/PrewittDiffCutOff/cutOff_" + to_string(CUTOFF) + "-filterSize_" + to_string(FILTER_SIZE) + ".txt");
	output_edge_co.open("../test_data/EdgeDetectionCutOff/cutOff_" + to_string(CUTOFF) + "-Depth_" + to_string(DEPTH) + ".txt");
	
	if(argc != __ARG_NUM__)
	{
		usage();
		return 0;
	}

	BitmapRawConverter inputFile(argv[1]);
	BitmapRawConverter outputFileSerialPrewitt(argv[1]);
	BitmapRawConverter outputFileParallelPrewitt(argv[1]);
	BitmapRawConverter outputFilePrewittParallelFor(argv[1]);
	BitmapRawConverter outputFileSerialEdge(argv[1]);
	BitmapRawConverter outputFileParallelEdge(argv[1]);

	unsigned int width, height;

	int test;
	
	width = inputFile.getWidth();
	height = inputFile.getHeight();

	int* outBufferSerialPrewitt = new int[width * height];
	int* outBufferParallelPrewitt = new int[width * height];
	int* outBufferPrewittParallelFor = new int[width * height];

	memset(outBufferSerialPrewitt, 0x0, width * height * sizeof(int));
	memset(outBufferParallelPrewitt, 0x0, width * height * sizeof(int));
	memset(outBufferPrewittParallelFor, 0x0, width * height * sizeof(int));

	int* outBufferSerialEdge = new int[width * height];
	int* outBufferParallelEdge = new int[width * height];

	memset(outBufferSerialEdge, 0x0, width * height * sizeof(int));
	memset(outBufferParallelEdge, 0x0, width * height * sizeof(int));

	// serial version Prewitt
	run_test_nr(1, &outputFileSerialPrewitt, argv[2], outBufferSerialPrewitt, width, height, output_prewitt_fs, output_prewitt_co);

	// parallel version Prewitt
	run_test_nr(2, &outputFileParallelPrewitt, argv[3], outBufferParallelPrewitt, width, height, output_prewitt_fs, output_prewitt_co);

	thresholdCutOff(outputFileSerialEdge.getBuffer(), width, height);

	// serial version special
	run_test_nr(3, &outputFileSerialEdge, argv[4], outBufferSerialEdge, width, height, output_edge_co , output_edge_co);

	thresholdCutOff(outputFileParallelEdge.getBuffer(), width, height);

	// parallel version special
	run_test_nr(4, &outputFileParallelEdge, argv[5], outBufferParallelEdge, width, height, output_edge_co, output_edge_co);

	// parallel version Prewitt using parallel_for
	run_test_nr(5, &outputFileParallelEdge, "outputPrewittParallelFor.bmp", outBufferPrewittParallelFor, width, height, output_prewitt_fs, output_prewitt_co);

	output_prewitt_fs.close();
	output_prewitt_co.close();
	output_edge_co.close();

	// verification
	cout << "Verification: ";
	test = memcmp(outBufferSerialPrewitt, outBufferParallelPrewitt, width * height * sizeof(int));

	if(test != 0)
	{
		cout << "Prewitt FAIL!" << endl;
	}
	else
	{
		cout << "Prewitt PASS." << endl;
	}

	test = memcmp(outBufferSerialEdge, outBufferParallelEdge, width * height * sizeof(int));

	if(test != 0)
	{
		cout << "Edge detection FAIL!" << endl;
	}
	else
	{
		cout << "Edge detection PASS." << endl;
	}

	// clean up
	delete outBufferSerialPrewitt;
	delete outBufferParallelPrewitt;

	delete outBufferSerialEdge;
	delete outBufferParallelEdge;

	return 0;
} 