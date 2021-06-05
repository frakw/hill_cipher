#include <iostream>
#include <string>
#include <math.h>

template<typename T>
class Matrix {
public:
	Matrix(int _row, int _col) {
		new_mem(_row, _col);
	}
	Matrix(int _row, int _col, T* input) {
		new_mem(_row, _col);
		load_data(input);
	}
	Matrix(const Matrix<T>& in) {//copy constructor
		new_mem(in.row, in.col);
		load_data(in.data);
	}
	~Matrix() {
		del_mem();
	}
	void print() {
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				std::cout << matrix[i][j] << '\t';
			}
			std::cout << std::endl;
		}
	}
	T* operator[](const int index)const {
		return matrix[index];
	}
	Matrix<T>& operator=(const Matrix<T>& input) {
		del_mem();
		new_mem(input.row, input.col);
		load_data(input.data);
		return *this;
	}
	const Matrix<T> operator*(T input) {
		Matrix<T> result(row, col);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				result[i][j] = this->matrix[i][j] * input;
			}
		}
		return result;
	}

	const Matrix<T> operator%(T input) {
		Matrix<T> result(row, col);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				result[i][j] = this->matrix[i][j] % input;
			}
		}
		return result;
	}
	const Matrix<T> operator*(const Matrix<T>& input) {
		if (this->col != input.row) {
			T null = -1;
			return Matrix<T>(1, 1, &null);
		}
		Matrix<T> result(this->row, input.col);
		for (int i = 0; i < this->row; i++) {
			for (int j = 0; j < input.col; j++) {
				result[i][j] = 0;
				for (int k = 0; k < this->col; k++) {
					result[i][j] += this->matrix[i][k] * input[k][j];
				}
			}
		}
		return result;
	}
	const Matrix<T> transpose() {
		Matrix<T> result(col, row);
		for (int i = 0; i < result.row; i++) {
			for (int j = 0; j < result.col; j++) {
				result[i][j] = matrix[j][i];
			}
		}
		return result;
	}
	Matrix<T> inverse() {
		if (this->col != this->row) {
			T null = -1;
			return Matrix<T>(1, 1, &null);
		}
		int n = this->row;
		Matrix<T> tmp(n - 1, n - 1);
		Matrix<T> result(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				sub_matrix(*this, tmp, i, j);
				int odd_count = (i % 2) + (j % 2);
				result[i][j] = (odd_count == 1 ? -1 : 1) * tmp.determinant();
			}
		}
		result = result.transpose();
		return result * (1.0 / (float)this->determinant());
	}
	Matrix<T> mod_inverse(T mod_num) {
		if (this->col != this->row) {
			T null = -1;
			return Matrix<T>(1, 1, &null);
		}
		int n = this->row;
		Matrix<T> tmp(n - 1, n - 1);
		Matrix<T> result(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				sub_matrix(*this, tmp, i, j);
				int odd_count = (i % 2) + (j % 2);
				result[i][j] = (odd_count == 1 ? -1 : 1) * tmp.determinant();
			}
		}
		result = result.transpose();
		result = result * extended_euclidean(this->determinant(), mod_num);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				result[i][j] = mod(result[i][j], mod_num);
			}
		}
		return result;
	}
	static void sub_matrix(Matrix<T>& src, Matrix<T>& dst, int row_num, int col_num) {
		if ((dst.row != (src.row - 1)) || (dst.col != (src.col - 1))) {
			dst.del_mem();
			dst.new_mem(src.row - 1, src.col - 1);
		}
		int dst_row_index = 0, dst_col_index = 0;
		for (int i = 0; i < src.row; i++) {
			for (int j = 0; j < src.col; j++) {
				if (i != row_num && j != col_num) {
					dst[dst_row_index][dst_col_index] = src[i][j];
					dst_col_index++;
					if (dst_col_index == dst.col) {
						dst_col_index = 0;
						dst_row_index++;
					}
				}
			}
		}
	}
	T determinant() {
		if (row != col) return -1;
		return _determinant(row);
	}
	int row, col;
private:
	void new_mem(int _row, int _col) {
		row = _row;
		col = _col;
		data = new T[row * col];
		matrix = new T * [row];
		for (int i = 0; i < row; i++) {
			matrix[i] = &data[i * col];
		}
	}
	void del_mem() {
		delete[] data;
		delete[] matrix;
	}
	void load_data(T* input) {
		for (int i = 0; i < row * col; i++) {
			data[i] = input[i];
		}
	}
	T _determinant(int n) {
		T result = 0;
		if (n == 1) {
			return matrix[0][0];
		}
		if (n == 2) {
			return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
		}
		Matrix<T> tmp(n - 1, n - 1);
		int sign = 1;
		for (int i = 0; i < n; i++) {
			sub_matrix(*this, tmp, 0, i);
			result += sign * matrix[0][i] * tmp.determinant();
			sign = -sign;
		}
		return result;
	}

	T mod(T x, T mod_num) {
		x = (int)x % (int)mod_num;
		if (x >= 0) return x;
		else return x + mod_num;
	}

	T extended_euclidean(T determinant, T mod_num) {
		determinant = mod(determinant, mod_num);
		for (T i = 2; i < mod_num; i++) {
			if ((int)(determinant * i) % (int)mod_num == 1) return i;
		}
		return -1;
	}

	T** matrix = nullptr;
	T* data = nullptr;
};

template <>
const Matrix<float> Matrix<float>::operator%(float input) {
	Matrix<float> result(row, col);
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i][j] = fmod(this->matrix[i][j], input);
		}
	}
	return result;
}

template<typename T>
Matrix<T> string_to_matrix(std::string input, int row, int col) {
	Matrix<T> result(row, col);
	int index = 0;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i][j] = input[index++] - 'a';
		}
	}
	return result;
}

template<typename T>
std::string matrix_to_string(Matrix<T> input) {
	std::string result;
	result.resize(input.row * input.col);
	int index = 0;
	for (int i = 0; i < input.row; i++) {
		for (int j = 0; j < input.col; j++) {
			result[index++] = input[i][j] + 'a';
		}
	}
	return result;
}


template<typename T>
std::string encrypt(std::string target, Matrix<T> key) {
	std::string result;
	if (target.length() % key.row != 0) {
		return "invalid message";
	}
	for (int i = 0; i < target.length(); i += 3) {
		Matrix<T> encrypt_matrix = key * string_to_matrix<T>(target.substr(i, 3), 3, 1);
		encrypt_matrix = encrypt_matrix % 26;
		result += matrix_to_string<T>(encrypt_matrix);
	}
	return result;
}


template<typename T>
std::string decrypt(std::string target, Matrix<T> key) {
	std::string result;
	if (target.length() % key.row != 0) {
		return "invalid message";
	}
	for (int i = 0; i < target.length(); i += 3) {
		Matrix<T> decrypt_matrix = key.mod_inverse(26) * string_to_matrix<T>(target.substr(i, 3), 3, 1);
		decrypt_matrix = decrypt_matrix % 26;
		result += matrix_to_string<T>(decrypt_matrix);
	}
	return result;
}

int main() {
	int key[9] =
	{
		17, 17,	5,
		21, 18, 21,
		2,	2,	19
	};
	Matrix<int> key_matrix(3, 3, key);
	std::string message = "paymoremoney";
	std::cout << "Original message : " << message << "\n\n--------------------------------\n\n";
	std::cout << "Key matrix : \n";
	key_matrix.print();
	std::string encrypted_message = encrypt(message, key_matrix);
	std::cout << "\nEncrypted message : " << encrypted_message << "\n\n--------------------------------\n\n";
	std::cout << "Inverse key matrix : \n";
	key_matrix.mod_inverse(26).print();
	std::string decrypted_message = decrypt(encrypted_message, key_matrix);
	std::cout << "\nDecrypted message : " << decrypted_message << std::endl;
	return 0;
}