#include "Methods.h"
#include <random>
#include <iostream>

const int RIGHT_BORDER = 100;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<int> distribution(0, RIGHT_BORDER);

SLE::SLE(int size = 0, double eps = 0.00001) : _size(size), _eps(eps)
{
	_matrix.resize(_size);
	_result.resize(METHODS_NUMBER);
	for (int i = 0; i < _size; ++i)
		_matrix[i].resize(_size + 1);
	for (int i = 0; i < METHODS_NUMBER; ++i)
		_result[i].resize(_size);
}

SLE::~SLE()
{
}

void SLE::resize(int newSize)
{
    _size = newSize;
    _matrix.resize(_size);
    for (int i = 0; i < _size; ++i)
    {
        _matrix[i].resize(_size + 1);
    }
    for (int i = 0; i < METHODS_NUMBER; ++i)
       _result[i].resize(_size);
}

//init
void SLE::ManualInput(vector<vector<double>> arr)
{
	for (int i = 0; i < _size; ++i)
	{
		for (int j = 0; j < _size + 1; ++j)
		{
			_matrix[i][j] = arr[i][j];
		}
	}
}

void SLE::AutomaticGeneration()
{
	for (int i = 0; i < _size; ++i)
	{
		for (int j = 0; j < _size + 1; ++j)
		{
			_matrix[i][j] = distribution(gen);
		}
	}
}

void SLE::GenerateSuitableMatrix()
{
	AutomaticGeneration();
	if (_symmetric)
	for (int i = 0; i < _size; ++i)
        for (int j = 0; j < _size; ++j)
			if (i > j)
				_matrix[j][i] = _matrix[i][j];
    if (_prevalence)
    for (int i = 0; i < _size; ++i)
    {
        double sum = 0;
        for (int j = 0; j < _size; ++j)
            sum += _matrix[i][j];
        _matrix[i][i] += sum;
    }
}

//methods
int SLE::Gauss()
{
	vector<vector<double>> arr(_size);
	for (int i = 0; i < _size; ++i)
	{
		arr[i].resize(_size + 1);
		for (int j = 0; j < _size + 1; ++j)
			arr[i][j] = _matrix[i][j];
	}
	steady_clock::time_point start = steady_clock::now();

	for(int col = 0; col < _size; ++col)
	{
		int maxIndex = 0;
		double maxValue = 0;

		for(int str = col; str < _size; ++str)
		{
			if (abs(arr[str][col]) > maxValue)
			{
				maxValue = abs(arr[str][col]);
				maxIndex = str;
			}
		}

		if (maxValue < DBL_EPSILON)
		{
			throw("solution with a given accuracy is impossible");
		}
		if (maxIndex != col)
		{
			std::swap(arr[col], arr[maxIndex]);
		}

		double tmp = arr[col][col];
		for (int i = col; i < _size + 1; ++i)
			arr[col][i] /= tmp;

		for (int i = 0 ; i < _size; ++i)
		{
			tmp = arr[i][col];
			if (abs(arr[i][col]) < _eps) continue;
			if (i == col) continue;
			for (int j = col; j < _size + 1; ++j)
				arr[i][j] -= arr[col][j] * tmp;
		}
	}
	for (int i = 0; i < _size; ++i)
		_result[GAUSS][i] = arr[i][_size];
	steady_clock::time_point end = steady_clock::now();
    int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
		(end - start).count();
	return elapsed_seconds;
}

int SLE::Cramer()
{
	vector<vector<double>> matrix(_size);
	for (int i = 0; i < _size; ++i)
	{
		matrix[i].resize(_size);
		for (int j = 0; j < _size; ++j)
			matrix[i][j] = _matrix[i][j];
	}
	steady_clock::time_point start = steady_clock::now();
	double determ = GetDeterminant(matrix);
	vector<double> col_determ(_size);
	for (int i = 0; i < _size; ++i)
	{
		for (int j = 0; j < _size; ++j)
		{
			matrix[j][i] = _matrix[j][_size];
		}
		col_determ[i] = GetDeterminant(matrix);
		for (int j = 0; j < _size; ++j)
			matrix[j][i] = _matrix[j][i];
	}
	for (int i = 0; i < _size; ++i)
		_result[CRAMER][i] = col_determ[i] / determ;
	steady_clock::time_point end = steady_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
		(end - start).count();
	return elapsed_seconds;
}

int SLE::FixedPointIteration()
{
	vector<double> betta(_size);
	vector<double> lastResult(_size);
	vector<double> currentResult(_size);
	for (int i = 0; i < _size; ++i)
		betta[i] = _matrix[i][_size] / _matrix[i][i];
	lastResult = betta;
	vector<vector<double>> matrix(_size);
	for (int i = 0; i < _size; ++i)
	{
		matrix[i].resize(_size);
		for (int j = 0; j < _size; ++j)
			matrix[i][j] = _matrix[i][j];

	}
	if (!checkDiagonalPredominance(matrix))
		throw("There is no diagonal preponderance");
	for (int i = 0; i < _size; ++i)
		for (int j = 0; j < _size; ++j)
			matrix[i][j] /= (-1) * _matrix[i][i];
	for (int i = 0; i < _size; ++i)
		matrix[i][i] = 0;

	steady_clock::time_point start = steady_clock::now();

	while (true)
	{
		double norm = 0;
		double normV = 0;
		for (int i = 0; i < _size; ++i)
		{
			double tmp = 0;
			for (int j = 0; j < _size; ++j)
				tmp += lastResult[j] * matrix[i][j];
			currentResult[i] = tmp + betta[i];
		}
		for (int i = 0; i < _size; ++i)
			for (int j = 0; j < _size; ++j)
				if (abs(matrix[i][j]) > norm)
					norm = abs(matrix[i][j]);
		for (int i = 0; i < _size; ++i)
			if (abs(lastResult[i] - currentResult[i]) > normV)
				normV = abs(lastResult[i] - currentResult[i]);
        /* if (normV <= (1 - norm) / norm * _eps )*/
        if (normV <= _eps)
			break;
		std::swap(lastResult, currentResult);
	}
	std::swap(currentResult, _result[FIXEDPOINT]);

	steady_clock::time_point end = steady_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
		(end - start).count();
	return elapsed_seconds;
}

int SLE::Relaxation(double omega)
{
	if (!isSymmetric(_matrix))
		throw("Matrix is not symmetric");
	if (!isPositiv(_matrix))
		throw("Matrix is not positiv");

	vector<double> lastResult(_size);
	for (int i = 0; i < _size; ++i)
		lastResult[i] = _matrix[i][_size];

	vector<double> currentResult(_size);
	steady_clock::time_point start = steady_clock::now();
	while (true)
	{
		double max = 0;
		for (int i = 0; i < _size; ++i)
		{
			double sum = 0;
			for (int j = 0; j < i; ++j)
				sum += (_matrix[i][j] * currentResult[j]);
			for (int j = i + 1; j < _size; ++j)
				sum += (_matrix[i][j] * lastResult[j]);
			currentResult[i] = ( (-1) * omega * sum + (1 - omega) * _matrix[i][i] * lastResult[i] + omega * _matrix[i][_size]) / _matrix[i][i];
		}
		for (int i = 0; i < _size; ++i)
			if (abs(currentResult[i] - lastResult[i]) > max)
				max = abs(currentResult[i] - lastResult[i]);
		if (max < _eps)
			break;
		std::swap(lastResult, currentResult);
	}
	if (omega == 1)
		std::swap(currentResult, _result[SEIDEL]);
	else
		std::swap(currentResult, _result[RELAXATION]);
	steady_clock::time_point end = steady_clock::now();
	int elapsed_seconds = std::chrono::duration_cast<std::chrono::microseconds>
		(end - start).count();
	return elapsed_seconds;
}

int SLE::Seidel()
{
	return Relaxation(1);
}



double SLE::GetDeterminant(vector<vector<double>> matrix)
{
	double determinant = 1;

    for (int col = 0; col < matrix.size(); ++col)
	{
		int pivot_index = -1;
		double pivot_value = 0;
        for (int j = col; j <  matrix.size(); ++j)
		{
			if (abs(matrix[j][col]) > pivot_value)
			{
				pivot_index = j;
				pivot_value = abs(matrix[j][col]);
			}
		}

		if (pivot_value < DBL_EPSILON)
		{
			return 0;
		}

		if (pivot_index != col)
		{
			std::swap(matrix[pivot_index], matrix[col]);
			determinant *= -1;
		}

        for (int j = col + 1; j <  matrix.size(); j++)
		{
			if (abs(matrix[j][col]) < DBL_EPSILON)
				continue;
			double tmp = matrix[j][col] / matrix[col][col];
            for (int k = col; k <  matrix.size(); k++)
			{
				matrix[j][k] -= matrix[col][k] * tmp;
			}
		}
		determinant *= matrix[col][col];
	}
	return determinant;
}

bool SLE::isPositiv(vector<vector<double>> matrix)
{
	vector<vector<double>> tmp;
	for (int i = 1; i <= matrix.size(); ++i)
	{
		tmp.resize(i);
		for (int j = 0; j < i; ++j)
		{
			tmp[j].resize(i);
			for (int k = 0; k < i; ++k)
				tmp[j][k] = matrix[j][k];
		}
        if (GetDeterminant(tmp) <= 0)
			return false;
	}
	return true;
}

bool SLE::isSymmetric(vector<vector<double>> matrix)
{
	for (int i = 0; i < matrix.size(); ++i)
		for (int j = 0; j < matrix.size(); ++j)
			if (matrix[i][j] != matrix[j][i])
				return false;
	return true;
}

bool SLE::checkDiagonalPredominance(vector<vector<double>> matrix)
{
	for (int i = 0; i < matrix.size(); ++i)
	{
        double sum = 0;
		for (int j = 0; j < matrix[i].size(); ++j)
		{
			sum += abs(matrix[i][j]);
		}
        if (sum > 2 * abs(matrix[i][i]))
			return false;
	}
	return true;
}

void SLE::resetResult()
{
    for (int i = 0; i < METHODS_NUMBER; ++i)
        for (int j = 0; j < _size; ++j)
            _result[i][j] = 0;
}
