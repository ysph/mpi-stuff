#include <omp.h>
#include <iostream>
#include <bitset>
#include <string>
#include <vector>
#include <time.h>
#include <cmath>
#include <cfenv>
using namespace std;

#define POWER_OF_TWO 3 // 2 ^ 3
#define HIGHEST_VALUE_IN_BYTE 256 // 11111111
#define BIGGEST_SINGLE_BIT 128 // 10000000
#define INFINITE 10000000

#pragma STDC FENV_ACCESS ON
double roundFloat(double x) {
    std::fenv_t save_env;
    std::feholdexcept(&save_env);
    double result = std::rint(x);
    if (std::fetestexcept(FE_INEXACT)) {
        auto const save_round = std::fegetround();
        std::fesetround(FE_TOWARDZERO);
        result = std::rint(std::copysign(0.5 + std::fabs(x), x));
        std::fesetround(save_round);
    }
    std::feupdateenv(&save_env);
    return result;
}

typedef unsigned int IMAGE; // bit array

class BitImage {
private:
	struct Mask {
		int A, B, C;
	};
	
	IMAGE *binimg;
	int height, width, size, pixels;
	int amount;
	
	int bit = 1, bitMask = 1;
	int byte = 0, bitInByte = 0;
	int components = 0;
	int numberOfThreads = 2;
	int** convex_hull = NULL;
	
	string matrix = "";
	std::vector<int*> result;
	
public:	
	BitImage(int x, int y) {
		if (x <= 0 || y <= 0) {
			std::cout << "Either size cannot be less or equal 0" << std::endl;
			exit(0);
		}	
		height = x;
		width = y;
		pixels = x * y;
		
		size = (x * y) >> POWER_OF_TWO;
		
		if (((x * y) % 8) != 0) {
			size += 1;
		}
		
		binimg = new IMAGE[size];
		
		for (int i = 0; i < size; i++) {
			binimg[i] = 0;
			std::bitset<8> bits(0);
			matrix += bits.to_string();
		}
		
		int redundant = (size << POWER_OF_TWO) - pixels;
		
		matrix.erase(matrix.size() - redundant);
	};

	void setUpAccess(int x, int y) {
		bit = ((x - 1) * width) + y;
		byte = 0;
		bitMask = 1;
		bitInByte = bit % 8;
		
		byte = (bit >> POWER_OF_TWO);

		if (bitInByte != 0) {
			bitMask = BIGGEST_SINGLE_BIT >> (bitInByte - 1);
			byte += 1;
		}
	}
	
	void setRandom() {
		srand(time(NULL)); // use current time
		
		int digit;
		matrix = "";
		for (int i = 0; i < size; i++) {
			digit = rand() % HIGHEST_VALUE_IN_BYTE;
			binimg[i] = digit;
			std::bitset<8> bits(digit);
			matrix += bits.to_string();
		}
	}
	
	void printAll() {
		printf("\nOur image");
		for (int i = 0; i < pixels; i++) {
			if ( (i % width) == 0 ) { putchar('\n'); }
			printf("%c", matrix[i]);
		}
		putchar('\n');
	}
	
	void setPixel(int x, int y) {
		if (x <= 0 || y <= 0 || pixels < (x * y)) {
			cout << "You are accessing " << (x * y) << "th bit from only " << pixels << " available, ";
			cout << "therefore it's a nonexistent bit, abort." << endl;
			return;
		}
		
		setUpAccess(x, y);
		
		cout << "Accessing " << byte << "th byte: " << std::bitset<8>(binimg[byte - 1]) << ", ";
		cout << "setting up " << bit << "th bit:\t" << std::bitset<8>(bitMask) << endl;
		
		binimg[byte - 1] |= bitMask;
		matrix[bit - 1] = '1';
	}
	
	void fillImage() {
		string s;
		matrix = "";
		for (int i = 0; i < size; i++) {
			binimg[i] = HIGHEST_VALUE_IN_BYTE - 1; // all bits set to 1's.
			auto s = to_string(binimg[i]);
			matrix += s;
		}
	}
	
	void fillHalf() {
		string s;
		matrix = "";
		for (int i = 0; i < size; i++) {
			if (!(i % 2)) continue;
			binimg[i] = HIGHEST_VALUE_IN_BYTE - 1; // all bits set to 1's.
			auto s = to_string(binimg[i]);
			matrix += s;
		}
	}
	
	bool getPixel(int x, int y) {
		setUpAccess(x, y);
		
		int nthByte = 7 - ((bit - 1) % 8);
		std::bitset<8> ourByte(binimg[byte - 1]);

		if (ourByte[nthByte]) return true;
		
		return false;
	}
	
	double turnPoint(int* p1, int* p2, int* p3) {
		double ax = p2[0] - p1[0];
  		double ay = p2[1] - p1[1];
  		double bx = p3[0] - p1[0];
  		double by = p3[1] - p1[1];
  		
  		return (ax * bx + ay * by) / (sqrt(ax * ax + ay * ay) * sqrt(bx * bx + by * by));
	}
	
	double betweenPoints(int* p1, int* p2) {
  		return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]));
	}
	void findConnectedComponents() {
		int kWidth = 0, kHeight = 0;
		string s; // convertion
		
		amount = 1;
		Mask localMask;
			
		for (int i = 1; i <= width; i++) {
			for (int j = 1; j <= height; j++) {
				kWidth = j - 1;
				if (kWidth <= 0) {
					kWidth = 1;
					localMask.B = 0;
				} else {
					localMask.B = getPixel(i, kWidth);
				}
				
				kHeight = i - 1;
      			if (kHeight <= 0) {
        			kHeight = 1;
        			localMask.C = 0;
      			} else {
        			localMask.C = getPixel(kHeight, j);
      			}
      			
      			localMask.A = getPixel(i, j);
      			
      			if (localMask.A == 0) {
      			} else if (localMask.B == 0 && localMask.C == 0) {
      				components += 1;
      				amount += 1;
          			s = to_string(amount);
          			matrix[((i - 1) * width) + (j - 1)] = s[0];	
        		} else if (localMask.B != 0 && localMask.C == 0) {
        			s = to_string(localMask.B);
          			matrix[((i - 1) * width) + (j - 1)] = s[0];
				} else if (localMask.B == 0 && localMask.C != 0) {
        			s = to_string(localMask.C);
          			matrix[((i - 1) * width) + (j - 1)] = s[0];
          		} else if (localMask.B != 0 && localMask.C != 0) {
        			s = to_string(localMask.B);
          			matrix[((i - 1) * width) + (j - 1)] = s[0];
          			
          			// find pixels labeled as C, then relabel them as B.
          			if (localMask.B != localMask.C) {
          				for (int a = 0; a < i; a++) {
              				for (int b = 0; b < height; b++) {
                				if (matrix[(a * width) + b] == localMask.C) {
                					s = to_string(localMask.B);
          							matrix[(a * width) + b] = s[0];	
                  				}
                  			}
                  		}
          			}
          			// fi
        		}
			}
		}
	}
	
	void makeConvex() {
		int tempo = 0;
		convex_hull = new int*[components];
		
		for (int i = 0; i < components; i++) {
			convex_hull[i] = new int[3];
		}

		for (int i = 0; i < width; i++) {
    		for (int j = 0; j < height; j++) {
    		
    			char access = matrix[(i * width) + j];
    			if (access != '0') {
    				convex_hull[tempo][0] = j;
    				convex_hull[tempo][1] = i;
    				
    				auto s = (int)access;
        			convex_hull[tempo][2] = s;
        			
    				if (tempo >= components - 1) break;
    				tempo += 1;
    			}
    		}
    	}
	}
	
	void Jarvis() {
		int m = 0, minind = 0;
		double mincos, cosine;
    	double len = 0, maxlen = 0;
    
		if (components == 1) {
    		result.push_back(convex_hull[0]);
    		return;
  		}
  		if (components == 2) {
    		result.push_back(convex_hull[0]);
    		result.push_back(convex_hull[1]);
    		return;
  		}
  		
  		double* first_elements = new double[2];
    	first_elements[0] = convex_hull[0][0];
    	first_elements[1] = convex_hull[0][1];
    	
    	for (int i = 1; i < components; i++) {
      		if (convex_hull[i][1] < first_elements[1]) {
        		m = i;
      		} else {
        		if ((convex_hull[i][1] == first_elements[1]) && (convex_hull[i][0] < first_elements[0])) {
          			m = i;
        		}
      		}
    	}
    	
    	result.push_back(convex_hull[m]);
    	
    	int* last = new int[2];
    	int* beforelast = new int[2];
    	
    	last = convex_hull[m];
    	beforelast[0] = convex_hull[m][0] - 2;
    	beforelast[1] = convex_hull[m][1];
    	
    	for(;;) {
    		mincos = 2;
      		for (int i = 0; i < components; i++) {
        		cosine = roundFloat(turnPoint(last, beforelast, convex_hull[i]) * INFINITE) / INFINITE;
        		if (cosine < mincos) {
          			minind = i;
          			mincos = cosine;
          			maxlen = betweenPoints(last, convex_hull[i]);
        		} else if (cosine == mincos) {
        			len = betweenPoints(last, convex_hull[i]);
          			if (len > maxlen) {
            			minind = i;
            			maxlen = len;
          			}
        		}
      		}
		
      		beforelast = last;
      		last = convex_hull[minind];
      		if (last == convex_hull[m])
        		break;
      		result.push_back(convex_hull[minind]);
    	}
	}
	
	void printResult() {
		int* temp = NULL;
		
		for (int i = 0; i < result.size(); i++) {
			temp = result[i];
		}
		putchar('\n');
		for (int i = 0; i < sizeof(temp)/sizeof(temp[0]); i++)
			cout << "(" << temp[i] << ")";
			
		putchar('\n');
	}
	
	void findConnectedComponents__OMP() {
		int kWidth = 0, kHeight = 0;
		string s; // converion
		
		amount = 1;
		Mask localMask;
			
		for (int i = 1; i <= width; i++) {
			for (int j = 1; j <= height; j++) {
				kWidth = j - 1;
				if (kWidth <= 0) {
					kWidth = 1;
					localMask.B = 0;
				} else {
					localMask.B = getPixel(i, kWidth);
				}
				
				kHeight = i - 1;
      			if (kHeight <= 0) {
        			kHeight = 1;
        			localMask.C = 0;
      			} else {
        			localMask.C = getPixel(kHeight, j);
      			}
      			
      			localMask.A = getPixel(i, j);
      			
      			if (localMask.A == 0) {
      			} else if (localMask.B == 0 && localMask.C == 0) {
      				components += 1;
      				amount += 1;
          			s = to_string(amount);
          			matrix[((i - 1) * width) + (j - 1)] = s[0];	
        		} else if (localMask.B != 0 && localMask.C == 0) {
          			matrix[((i - 1) * width) + (j - 1)] = (char)localMask.B;
				} else if (localMask.B == 0 && localMask.C != 0) {
        			s = (char)localMask.C;
          			matrix[((i - 1) * width) + (j - 1)] = s[0];
          		} else if (localMask.B != 0 && localMask.C != 0) {
          			matrix[((i - 1) * width) + (j - 1)] = (char)localMask.B;
          			
          			if (localMask.B != localMask.C) {
          			// find pixels labeled as C, then relabel them as B.
          			int a;
          			
					#pragma omp parallel private(a) num_threads(numberOfThreads) 
						{
					#pragma omp for
							for (a = 0; a < i; a++) {
              					for (int b = 0; b < height; b++) {
              						#pragma omp critical
        							{
                						if (matrix[(a * width) + b] == localMask.C) {
          									matrix[(a * width) + b] = (char)localMask.B;	
                  						}
                  					}
                  				}
                  			}
          				}
          			}
          			// fi
        		}
			}
		}
	}
	
	void Jarvis__OMP() {
		int id, delta, rest; // for omp
		int m = 0, minind = 0;
		double mincos, cosine;
    	double len = 0, maxlen = 0;
    	int* last;
    	int* beforelast;
    	
		std::vector<std::vector<int*>> local_result(numberOfThreads);
      		
		if (components == 1) {
    		result.push_back(convex_hull[0]);
    		return;
  		}
  		if (components == 2) {
    		result.push_back(convex_hull[0]);
    		result.push_back(convex_hull[1]);
    		return;
  		}
  		
  		result.push_back(convex_hull[m]);
  		
    	for (int i = 0; i < numberOfThreads; i++)
      		local_result[i].push_back(convex_hull[m]);
      
  		// OMP	
  		#pragma omp parallel private(id, delta, rest, mincos, cosine, minind, maxlen, last, beforelast, len) num_threads(numberOfThreads) 
  		{
  			delta = components / numberOfThreads;
			id = omp_get_thread_num();
		
			maxlen = minind = rest = 0;

      		if (id == numberOfThreads - 1) {
        		rest = components % numberOfThreads;
      		}
			
      		last = new int[2];
      		beforelast = new int[2];
      		
      		last = convex_hull[m];
      		beforelast[0] = convex_hull[m][0] - 2;
      		beforelast[1] = convex_hull[m][1];
		
      		for(;;) {
        		mincos = 2;
        		for (int i = id * delta; i < id * delta + delta + rest; i++) {
          			cosine = roundFloat(turnPoint(last, beforelast, convex_hull[i]) * INFINITE) / INFINITE;
          			if (cosine < mincos) {
            			minind = i;
            			mincos = cosine;
            			maxlen = betweenPoints(last, convex_hull[i]);
          			} else if (cosine == mincos) {
          			
            			len = betweenPoints(last, convex_hull[i]);
            			if (len > maxlen) {
              				minind = i;
              				maxlen = len;
            			}
          			}
        		}
		
        		if (id != 0) {
          			cosine = roundFloat(turnPoint(last, beforelast, convex_hull[0]) * INFINITE) / INFINITE;
          			if (cosine < mincos) {
            			minind = 0;
            			mincos = cosine;
            			maxlen = betweenPoints(last, convex_hull[0]);
          			}
          			if (cosine == mincos) {
          			
            			len = betweenPoints(last, convex_hull[0]);
            			if (len > maxlen) {
              				minind = 0;
              				maxlen = len;
            			}
          			}
        		}
		
        		beforelast = last;
        		last = convex_hull[minind];
        		
        		if (last == convex_hull[m])
          			break;
        		local_result[id].push_back(convex_hull[minind]);
      		}
  		}
  		// END OF OMP
  		
  		std::vector<int*> finale_local;
    	for (int i = 0; i < numberOfThreads; i++) {
      		int finalSize = local_result[i].size();
      		
      		for (int j = 0; j < finalSize; j++)
        		finale_local.push_back(local_result[i][j]);
    		}
		
    		int s = finale_local.size();
		
    		last = new int[2];
    		beforelast = new int[2];
    		last = convex_hull[m];
    		beforelast[0] = convex_hull[m][0] - 2;
    		beforelast[1] = convex_hull[m][1];
    		
    		for(;;) {
      			mincos = 2;
      			for (int i = 0; i < s; i++) {
        			cosine = roundFloat(turnPoint(last, beforelast, finale_local[i]) * INFINITE) / INFINITE;
        			
        			if (cosine < mincos) {
          				minind = i;
          				mincos = cosine;
          				maxlen = betweenPoints(last, finale_local[i]);
        			}
        			
        			if (cosine == mincos) {
          				len = betweenPoints(last, finale_local[i]);
          				if (len > maxlen) {
            				minind = i;
            				maxlen = len;
          				}
        			}
      			}
			
      			beforelast = last;
      			last = finale_local[minind];
      			if (last == finale_local[m]) break;
      			result.push_back(finale_local[minind]);
    		}
	}
	
	~BitImage() {
		matrix.clear();
	}
};

int main() {
	int height, width;
	
	cout << "We represent our binary image as a packed one-dimensional array in which each pixel is a bit." << endl;
	/* user input */
	/*cout << "Enter height, then width: ";
	cin >> height >> width;
	putchar('\n');*/
	
	/* predefined input */
	width = 4500;
	height = 4500;

	/* create binary image (object) */
	BitImage be(height, width);
	BitImage ae(height, width);
	
	be.setRandom(); // set random image
	ae = be;
	
	//start measuring time
	clock_t begin = clock();
	//doing job
	be.findConnectedComponents__OMP();
	be.makeConvex();
	be.Jarvis__OMP();
	//stop measuring time
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//end.
	cout << "(OMP) Time spent on finding connected components and making the convex hull is: " << time_spent << " seconds" << endl;
	
	//start measuring time
	begin = clock();
	//doing job
	ae.findConnectedComponents();
	ae.makeConvex();
	ae.Jarvis();
	//stop measuring time
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	cout << "(Seq) Time spent on finding connected components and making the convex hull is: " << time_spent << " seconds" << endl;
	be.printResult();
	
	return 0;
}
