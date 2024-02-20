#include "s21_matrix_plus.h"

//base constructor witout arguments
S21Matrix::S21Matrix() {
            rows_ = 0;
            cols_ = 0;

            matrix_ = new double*[rows_ * cols_];
            for (int i = 0; i < rows_; i++)
                matrix_[i] = new double[cols_];

            for (int n = 0; n < this->rows_; n++)
                for (int m = 0; m < this->cols_; m++)
                    this->matrix_[n][m] = 0;

            //std::cout << "Call constructor default " << this << std::endl;
        }


//constructor with paranetrs
S21Matrix::S21Matrix(int rows, int cols) {
            rows_ = rows;
            cols_ = cols;

            matrix_ = new double*[rows_ * cols_];
            for (int i = 0; i < rows_; i++)
                matrix_[i] = new double[cols_];

            for (int n = 0; n < this->rows_; n++)
                for (int m = 0; m < this->cols_; m++)
                    this->matrix_[n][m] = 0;

            //std::cout << "Call constructor with parametrs " << this << std::endl;
        }

S21Matrix::S21Matrix(const S21Matrix& other) : S21Matrix(other.rows_, other.cols_) {
            for (int n = 0; n < this->rows_; n++)
                for (int m = 0; m < this->cols_; m++)
                    this->matrix_[n][m] = other.matrix_[n][m];
        }

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {

  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

//destructor
S21Matrix::~S21Matrix() {
            if (matrix_ ) {
                for (int i = 0; i < rows_; i++)
                    delete matrix_[i];
                delete []matrix_;
            }

            //std::cout << "Call destructor " << this << std::endl;
        }



S21Matrix S21Matrix::operator+(const S21Matrix &other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_  ) throw std::logic_error("Size of matrix mast be equal\n");        // сюда добавить исключение

    S21Matrix result(*this);
    result.SumMatrix(other);

    return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_  ) throw std::logic_error("Size of matrix mast be equal\n");        // сюда добавить исключение

    S21Matrix result(*this);
    result.SubMatrix(other);

    return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
    if (this->cols_ != other.rows_) throw std::logic_error("The number of columns of the first matrix is not equal to the number of rows of the second matrix.");

    S21Matrix result(this->rows_, other.cols_);
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        double number = 0;
        for (int n = 0; n < this->cols_; n++)
          number += this->matrix_[i][n] * other.matrix_[n][j];
        result.matrix_[i][j] = number;
      }
    }
    return result;
}

S21Matrix S21Matrix::operator*(const double num) {
    S21Matrix result(this->rows_, this->cols_);
    for (int i = 0; i < result.rows_; i++) {
      for (int j = 0; j < result.cols_; j++)
        result.matrix_[i][j] = this->matrix_[i][j] * num;
    }
    return result;
}

bool S21Matrix::operator==(const S21Matrix &other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_  ) return false;

    for (int i = 0; i < this->rows_; i++)
        for (int j = 0; j < this->cols_; j++)
            if (abs(this->matrix_[i][j] - other.matrix_[i][j]) > 0.001) return false;
        
    return true;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
    if (this->matrix_ ) {
                for (int i = 0; i < rows_; i++)
                    delete this->matrix_[i];
                delete []this->matrix_;
            }

    this->rows_ = other.rows_;
    this->cols_ = other.cols_;

    for (int i = 0; i < this->rows_; i++)
        for (int j = 0; j < this->cols_; j++)
            this->matrix_[i][j] = other.matrix_[i][j];

    return *this;
}


S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
    this->SumMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
    this->SubMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
    this->MulMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator*=(const double num) {
    this->MulNumber(num);
    return *this;
}

double &S21Matrix::operator()(int row, int col) & {
  return const_cast<double &>(this->matrix_[row][col]);
}


bool S21Matrix::EqMatrix(const S21Matrix& other) {
    return *this == other;
}




void S21Matrix::SumMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_  ) throw std::logic_error("Size of matrix mast be equal\n");        // сюда добавить исключение

    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] = this->matrix_[i][j] + other.matrix_[i][j];
    }

}

void S21Matrix::SubMatrix(const S21Matrix& other) {
    if (this->rows_ != other.rows_ || this->cols_ != other.cols_  ) throw std::logic_error("Size of matrix mast be equal\n");        // сюда добавить исключение

    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] = this->matrix_[i][j] - other.matrix_[i][j];
    }
}

void S21Matrix::MulNumber(const double num) {
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < this->cols_; j++)
        this->matrix_[i][j] = this->matrix_[i][j] * num;
    }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
    if (this->cols_ != other.rows_) throw std::logic_error("The number of columns of the first matrix is not equal to the number of rows of the second matrix.");

    S21Matrix result(this->rows_, other.get_cols());
    for (int i = 0; i < this->rows_; i++) {
      for (int j = 0; j < other.cols_; j++) {
        double number = 0;
        for (int n = 0; n < this->cols_; n++)
          number += this->matrix_[i][n] * other.matrix_[n][j];
        result.matrix_[i][j] = number;
      }
    }

    for (int i = 0; i < this->rows_; i++)
                    delete this->matrix_[i];
                delete []matrix_;

            this->rows_ = result.rows_;
            this->cols_ = result.cols_;

            this->matrix_ = new double*[result.rows_ * result.cols_];
            for (int i = 0; i < rows_; i++)
                matrix_[i] = new double[result.cols_];

            for (int n = 0; n < result.rows_; n++)
                for (int m = 0; m < result.cols_; m++)
                    this->matrix_[n][m] = result.matrix_[n][m];

}

S21Matrix S21Matrix::Transpose() {
    S21Matrix result(this->cols_, this->rows_);
    for (int i = 0; i < this->cols_; i++)
        for (int j = 0; j < this->rows_; j++) 
            result.matrix_[i][j] = this->matrix_[j][i];
    return result;
}

S21Matrix S21Matrix::CalcComplements() {
    if (this->rows_ != this->cols_) throw std::logic_error("The matrix is not square.");

    int size = this->rows_;
    S21Matrix result(size, size);

    for (int a = 0; a < size; a++) {
      for (int b = 0; b < size; b++) {
        //делаем вычеркнутую матрицу
        S21Matrix small_matrix(size-1, size-1);
        int i_small = 0;
        int j_small = 0;
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < size; j++) {
            if (i == a || j == b) continue;
            small_matrix.matrix_[i_small][j_small] = this->matrix_[i][j];
            j_small++;
            if (j_small == size - 1) {
              i_small++;
              j_small = 0;
            }
          }
        }
        //сделали

        //находим минор и алгебраическое дополнение
        double minor = 0;
        minor = small_matrix.Determinant();                              //this
        double algebraic_complement = pow(-1, (a + 1 + b + 1)) * minor;
        result.matrix_[a][b] = algebraic_complement;
      }
    }
    return result;
}

//возвращает определитель матрицы 2 на 2
int S21Matrix::s21_2x2_determinant(S21Matrix A) {
  int res =
      (A.matrix_[0][0] * A.matrix_[1][1]) - (A.matrix_[0][1] * A.matrix_[1][0]);
  return res;
}

double S21Matrix::Determinant() {
    if (this->rows_ != this->cols_) throw std::logic_error("The matrix is not square.");
    double result = 0;

  if (this->rows_ == 1 && this->cols_ == 1)
    result = this->matrix_[0][0];
  else if (this->rows_ == 2 && this->cols_ == 2)
    result = S21Matrix::s21_2x2_determinant(*this);
  else {
    //здесь нужно сделать рекурсивный поиск определителей, пока размер не станет
    // 2 на 2. считаем детернминант как сумму (элемент первой строки умноженный
    // на алгебраическое дополнение этого элемента) детерминант = сумма
    //элементов первой строки, умноженных на их алгебраические дополнения
    //детерминант = сумма элементов первой строки, умноженных на (минор *
    //(-1)^(i+j)) детерминант = сумма элементов первой строки, умноженных на
    // ((детерминант матрицы с вычеркнутой строкой i и вычеркнутым столбцом j) *
    // (-1)^(i+j))
    int size = this->rows_;

    //сумма
    for (int n = 0; n < size; n++) {
      //делаем вычеркнутую матрицу
      S21Matrix small_matrix(size -1, size -1);
      int i_small = 0;
      int j_small = 0;
      for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
          if (i == 0 || j == n) continue;
          small_matrix.matrix_[i_small][j_small] = this->matrix_[i][j];
          j_small++;
          if (j_small == size - 1) {
            i_small++;
            j_small = 0;
          }
        }
      }

      //находим минор и алгебраическое дополнение
      double minor = 0;
      minor = small_matrix.Determinant();
      double algebraic_complement = pow(-1, (1 + n + 1)) * minor;

      //прибавляем к результату число, равное произведению текущего элемента на
      //минор этого элемента
      result += this->matrix_[0][n] * algebraic_complement;

    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
    if (this->rows_ != this->cols_) throw std::logic_error("The matrix is not square.");                //exeption

    int size = this->rows_;
    S21Matrix result(size, size);

    if (size == 1) {
      if (this->matrix_[0][0] == 0) throw std::logic_error("Matrix determinant is 0.");              // exeption
      result.matrix_[0][0] = 1 / this->matrix_[0][0];
    } else {
    double determinant = 0;
    determinant = this->Determinant();

    if (determinant == 0)
      throw std::logic_error("Matrix determinant is 0.");                                         //exeption
    else {
      S21Matrix algebraic_complements = this->CalcComplements();
      S21Matrix transpose_algebraic_complements = algebraic_complements.Transpose();

      double coefficient = 1 / determinant;
      transpose_algebraic_complements.MulNumber(coefficient);
      return transpose_algebraic_complements;
    }
  }
  return result;
}

void S21Matrix::set_rows(int new_rows) {
  if (new_rows < 0) {
    throw std::length_error("matrix rows count must be non-negative");
  }

  if (new_rows != rows_) {

    S21Matrix tmp{new_rows, cols_};
    int min = std::min(rows_, new_rows);
    for (int i = 0; i < min; ++i) {
      for (int j = 0; j < cols_; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }

    for (int i = 0; i < this->rows_; i++)
      delete this->matrix_[i];
        delete []matrix_;

    this->matrix_ = new double*[new_rows * cols_];
            for (int i = 0; i < rows_; i++)
                matrix_[i] = new double[cols_];

            for (int n = 0; n < new_rows; n++)
                for (int m = 0; m < cols_; m++)
                    this->matrix_[n][m] = tmp.matrix_[n][m];
  }
}

void S21Matrix::set_cols(int new_cols) {
  if (new_cols < 0) {
    throw std::length_error("matrix columns count must be non-negative");
  }

  if (new_cols != cols_) {

    S21Matrix tmp{rows_, new_cols};
    int min = std::min(cols_, new_cols);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < min; ++j) {
        tmp(i, j) = (*this)(i, j);
      }
    }

    for (int i = 0; i < this->rows_; i++)
      delete this->matrix_[i];
        delete []matrix_;

    this->matrix_ = new double*[rows_ * new_cols];
            for (int i = 0; i < rows_; i++)
                matrix_[i] = new double[new_cols];

            for (int n = 0; n < rows_; n++)
                for (int m = 0; m < new_cols; m++)
                    this->matrix_[n][m] = tmp.matrix_[n][m];

  }
}


/*
int main() {

  S21Matrix mat_1(2, 2);

  mat_1.setValue(0, 0, 2);
  mat_1.setValue(0, 1, 2);

  mat_1.setValue(1, 0, 3);
  mat_1.setValue(1, 1, -2);

  S21Matrix mat_2(2, 2);

  mat_2.setValue(0, 0, 0.2);
  mat_2.setValue(0, 1, 0.2);

  mat_2.setValue(1, 0, 0.3);
  mat_2.setValue(1, 1, -0.2);

  S21Matrix buf(mat_1.InverseMatrix());

  mat_1.Print();
  mat_2.Print();

  buf.Print();

  bool result = mat_2.EqMatrix(buf);

  std::cout << "RES " << result << std::endl;


    return 0;
}

*/
















