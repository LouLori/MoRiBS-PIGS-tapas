#include <fstream>
#include <iostream>

using namespace std;

int main() {
	ifstream f("bb");

	int size_total, block_numb, size;
	f >> size_total>>block_numb;
	f >> size;

	int *vectors = new int[4*size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 4; j++) {
			f >> vectors[i*4+j];
		}
	}

	double *matrix = new double[size*size];
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			f >> matrix[i*size+j];
		}
	}

	f.close();

	cout << size << endl;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 4; j++) {
			cout << vectors[i*4+j] << " ";
		}
		cout << endl;
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << matrix[i*size+j] << " ";
		}
		cout << endl;
	}

	delete[] vectors;
	delete[] matrix;

	return 0;
}
