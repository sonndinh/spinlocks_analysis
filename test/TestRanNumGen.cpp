#include <chrono>
#include <random>
#include <iostream>

using namespace std;
#define MEAN 500
#define STDEV 200

void test_prng_time_based_seed();
void test_prng_fixed_seed();

int main(int argc, char** argv) {
	test_prng_time_based_seed();

	for (int i=0; i<5; i++) {
		test_prng_fixed_seed();
	}
	
	return 0;
}

void test_prng_time_based_seed() {	
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	cout << "Seed: " << seed << endl;
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
}

void test_prng_fixed_seed() {
	unsigned int seed = 5;
	minstd_rand generator(seed);
	//	knuth_b generator(seed);
	normal_distribution<double> distribution(MEAN, STDEV);

	for (int i=0; i<10; i++) {
		cout << distribution(generator) << " ";
	}
	cout << endl;
}
