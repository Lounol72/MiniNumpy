#Repertoires
BUILD_DIR = build
SRC_DIR = src
APP_DIR = apps
INCLUDE_DIR = include
TEST_DIR = tests
THIRDPARTY_DIR = third_party

# Noms des exécutables
TARGET = main
TARGET_VECTEUR = main_vector
TARGET_TEST = main_tests

# Fichiers source et objets
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
APPS = $(wildcard $(APP_DIR)/*.cpp)
DEPS = $(wildcard $(BUILD_DIR)/*.d)

# Tests source et objets
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -Werror -g -fsanitize=address,undefined -MMD -MP

EXES = $(TARGET) $(TARGET_VECTEUR) $(TARGET_TEST)
LIB = $(BUILD_DIR)/libmnp.a

all: $(EXES)

test: $(TARGET_TEST)
	./$(TARGET_TEST)

vecteur: $(TARGET_VECTEUR)
	./$(TARGET_VECTEUR)

$(BUILD_DIR):
	@mkdir -p $@

$(LIB): $(OBJECTS)
	ar -crs $@ $^

$(EXES): %: $(BUILD_DIR)/%.o $(LIB)
	@echo "Edition des lines de $@"
	@$(CXX) $(CXXFLAGS) $< -L $(BUILD_DIR) -lmnp -o $@

vpath %.cpp $(SRC_DIR) $(APP_DIR) $(TEST_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@echo "Compilation de $@"
	@$(CXX) $(CXXFLAGS) -I $(THIRDPARTY_DIR) -I $(INCLUDE_DIR) -c $< -o $@

clean:
		rm -rf $(BUILD_DIR) docs/carnet/build/*
		rm -f $(EXES)

doc:
	latexmk -pdf -cd -outdir=. -auxdir=build docs/carnet/carnet_maths.tex

.PHONY : all clean test vecteur

-include $(DEPS)
