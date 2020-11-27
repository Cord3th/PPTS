#include <ctime>
#include <iostream>
#include <fstream>

using namespace std;

int main() {
	ofstream output_file("my.map");
	int a[512];

	srand(time(0));
	for (int i = 0; i < 512; i++) {
		a[i] = i;
	}

	for (int i = 0 ; i < 512; i++) {
		int k = rand() % 512, temp = a[i];
		a[i] = a[k];
		a[k] = temp;
	}

	for (int i = 0; i < 512; i++) {
		output_file << (a[i] / 8 / 8 ) % 8 << ' '
				 	<< (a[i] / 8) % 8 << ' '
					<< a[i] % 8 << ' ' << 0 << endl;
	}

	output_file.close();

	return 0;
}
