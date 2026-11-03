CXX=g++
CXXFLAGS= -std=c++20 -Wall -Wextra -Werror -g -fsanitize=address,undefined
TARGET= main
TARGET_TEST= main_test
INCLUDES= include 
SOURCES_DIR=src
TEST_DIR=tests
THIRDPARTY_DIR= third_party

all : $(TARGET) $(TARGET_TEST)

test :  $(TARGET_TEST)
	./$(TARGET_TEST)

$(TARGET): main.o  version.o
		$(CXX) $(CXXFLAGS) $^ -o $@

$(TARGET_TEST): main_test.o  version.o vector.o
		$(CXX) $(CXXFLAGS) $^ -o $@

main.o: apps/hello.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLUDES) -c $< -o $@

main_test.o: tests/main_tests.cpp
	$(CXX) $(CXXFLAGS) -I $(THIRDPARTY_DIR) -I $(INCLUDES) -c $< -o $@

version.o: src/version.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLUDES) -c $< -o $@
vector.o: src/Vector.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLUDES) -c $< -o $@

clean:
		rm -f *.o $(TARGET) $(TARGET_TEST)

.PHONY : all clean test

