#include <iostream>
#include <cmath>

//Нужно доделать перегрузку оператора индексации, проверить работу оператора копирования
//Сделать функции из задания (все кроме трёх последних с помощью операторов)
//Сделать исключения
//Залить это всё в тесты
class S21Matrix {
    private:
        int rows_, cols_;
        double **matrix_;

    public:

        //Базовый конструктор, инициализирующий матрицу некоторой заранее заданной размерностью.
        S21Matrix();
        //Параметризированный конструктор с количеством строк и столбцов.
        S21Matrix(int rows, int cols);  
        //Конструктор копирования
        S21Matrix(const S21Matrix& other);
        //Конструктор переноса.    Вот тут вообще хуй знает, правильно ли
        S21Matrix(S21Matrix&& other);
        //Деструктор
        ~S21Matrix();

        void Print() {
            for (int n = 0; n < this->rows_; n++) {
                for (int m = 0; m < this->cols_; m++) {
                    std::cout << this->matrix_[n][m] << '\t';
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }

        void setValue(int a, int b, double n) {
            this->matrix_[a][b] = n;
        }

        int s21_2x2_determinant(S21Matrix A);

        //функции из задания
        bool EqMatrix(const S21Matrix& other);
        void SumMatrix(const S21Matrix& other);
        void SubMatrix(const S21Matrix& other);
        void MulNumber(const double num);
        void MulMatrix(const S21Matrix& other);
        S21Matrix Transpose();
        S21Matrix CalcComplements();
        double Determinant();
        S21Matrix InverseMatrix();


        int get_cols() const { return cols_; }
        int get_rows() const { return rows_; }
        void set_rows(int new_rows);
        void set_cols(int new_cols);

        //перегрузка операторов
        S21Matrix operator+(const S21Matrix& other);
        S21Matrix operator-(const S21Matrix& other);
        S21Matrix operator*(const S21Matrix& other);
        S21Matrix operator*(const double num);
        bool operator==(const S21Matrix& other);
        S21Matrix& operator=(const S21Matrix& other);
        S21Matrix& operator+=(const S21Matrix& other);
        S21Matrix& operator-=(const S21Matrix& other);
        S21Matrix& operator*=(const S21Matrix& other);
        S21Matrix& operator*=(const double num);
        double& operator()(int row, int col) &;
        //double& operator[][](int i, int j);

};