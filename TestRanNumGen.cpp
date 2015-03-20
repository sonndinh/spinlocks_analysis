#include <chrono>
#include <random>
#include <iostream>

using namespace std;
#define MEAN 500
#define STDEV 200

int main(int argc, char** argv) {
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	cout << "My seed is: " << seed << endl;
	knuth_b generator(seed);

	// provide mean and standard deviation
	normal_distribution<double> distribution(MEAN,STDEV);
	
	double tmp;
	int count = 0;
	for (int i=0; i<1000; i++) {
		tmp = distribution(generator);
		// eliminate those out of range [1,1000]
		if (tmp < 1 || tmp > 1000) {
			cout << "Out of range num: " << tmp << endl;
			count++;
		}
		//		cout << distribution(generator) << endl;
	}

	cout << "Number of out-of-range num is: " << count << endl;
	return 0;
}
