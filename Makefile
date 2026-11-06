#Repertoires
SRC_DIR= src
APP_DIR= apps
INCLUDE_DIR= include
TEST_DIR= tests
THIRDPARTY_DIR= third_party

# Noms des exécutables
TARGET= main
TARGET_VECTEUR= main_vector
TARGET_TEST= main_tests

# Fichiers source et objets
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SRCS:.cpp=.o)
APPS = $(wildcard $(APP_DIR)/*.cpp)

# Tests source et objets
CXX=g++
CXXFLAGS= -std=c++20 -Wall -Wextra -Werror -g -fsanitize=address,undefined

all : $(TARGET) $(TARGET_TEST)  $(TARGET_VECTEUR)

test :  $(TARGET_TEST)
	./$(TARGET_TEST)

vecteur : $(TARGET_VECTEUR)
	./$(TARGET_VECTEUR)

$(TARGET): main.o version.o
	$(CXX) $(CXXFLAGS)  $^ -o $@

$(TARGET_TEST): main_tests.o  version.o Vector.o
		$(CXX) $(CXXFLAGS) $^ -o $@

$(TARGET_VECTEUR): main_vector.o Vector.o
		$(CXX) $(CXXFLAGS) $^ -o $@

%.o: $(TEST_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS) -I $(THIRDPARTY_DIR) -I $(INCLUDE_DIR) -c $< -o $@

%.o: $(APP_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS) -I $(THIRDPARTY_DIR) -I $(INCLUDE_DIR) -c $< -o $@

%.o: $(SRC_DIR)/%.cpp
	@$(CXX) $(CXXFLAGS)  -I $(INCLUDE_DIR) -c $< -o $@

clean:
		rm -f *.o $(TARGET) $(TARGET_TEST) $(TARGET_VECTEUR)
		rm -f docs/carnet/build/*

doc:
	latexmk -pdf -cd -outdir=. -auxdir=build docs/carnet/carnet_maths.tex

.PHONY : all clean test vecteur

