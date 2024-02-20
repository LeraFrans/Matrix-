CC = g++ -std=c++17 -Wall -Werror -Wextra -Wpedantic
SOURCE = matrix.cpp
OBJECT = $(SOURCE:.cc=.o)
TEST_FLAGS =-lgtest

all: clean test

s21_matrix_oop.a:$(SOURCE)
	$(CC) -c $(SOURCE)
	@ar rcs s21_matrix_oop.a $(OBJECT)

test: clean unit_tests.cc s21_matrix_oop.a
	$(CC) unit_tests.cc matrix.cpp $(TEST_FLAGS) -o test
	./test



clang:
	clang-format -style=Google -i *.cc *.cpp  *.h

leaks_test: clean test
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	valgrind --leak-check=full --show-reachable=yes --show-leak-kinds=all --track-origins=yes --log-file=RESULT_VALGRIND.txt ./test
endif

clean:
	@rm -rf *.o *.a report.info *.gcda *.gcno report RESULT_VALGRIND.txt test