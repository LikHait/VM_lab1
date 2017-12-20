#pragma once
#include <chrono>
#include <vector>
using std::vector;

using std::chrono::steady_clock;
const int METHODS_NUMBER = 5;
const int GAUSS = 0;
const int CRAMER = 1;
const int FIXEDPOINT = 2;
const int SEIDEL = 3;
const int RELAXATION = 4;

class SLE
{
public:
	SLE(int size, double eps);
	~SLE();

    void resize(int newSize);

	//init
	void ManualInput(vector<vector<double>> arr);
	void AutomaticGeneration();
	void GenerateSuitableMatrix();
    void SetEpsilon(double eps) { _eps = eps; }

	//get
	vector<vector<double>> GetMatrix() const { return _matrix; }
	int GetSize() const { return _size; }
    vector<vector<double>> GetResult() const { return _result; }

	//methods
    int Gauss();
	int Cramer();
	int FixedPointIteration();
	int Relaxation(double omega = 1.5);
	int Seidel();

	double GetDeterminant(vector<vector<double>> matrix);
	bool isPositiv(vector<vector<double>> matrix);
	bool isSymmetric(vector<vector<double>> matrix);
	bool checkDiagonalPredominance(vector<vector<double>> matrix);
    void resetResult();

	//change bool
	void prevalenceToTrue() { _prevalence = true; }
	void symmetricToTrue() { _symmetric = true; }
    void prevalenceToFalse() { _prevalence = false; }
    void symmetricToFalse() { _symmetric = false; }
private:
	bool _prevalence = false;
	bool _symmetric = false;
	double _eps;
	int _size;
	vector<vector<double>> _matrix;
	vector<vector<double>> _result;
};
