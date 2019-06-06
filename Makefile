EXEC=test
SRC=.

FLAGS=-Wall -O3 -fopenmp
LINK=
INCLUDE=

BUILD_FOLDER=build2
# SRC_CPP=$(wildcard $(SRC)/*.cpp)
SRC_CPP=./test.cpp ./lines.cpp ./stb_image.cpp
OBJ_CPP=$(patsubst $(SRC)/%.cpp, $(BUILD_FOLDER)/%.o, $(SRC_CPP))

DEBUG_FOLDER=debug
DEB=$(DEBUG_FOLDER)/$(EXEC)
DEB_C=$(patsubst $(SRC)/%.c, $(DEBUG_FOLDER)/%.o, $(SRC_C))
DEB_CPP=$(patsubst $(SRC)/%.cpp, $(DEBUG_FOLDER)/%.o, $(SRC_CPP))

all: $(BUILD_FOLDER) $(EXEC)

$(BUILD_FOLDER):
	mkdir $(BUILD_FOLDER)

$(EXEC): $(OBJ_C) $(OBJ_CPP)
	g++ $(FLAGS) $^ -o $@ $(LINK)

$(BUILD_FOLDER)/%.o: $(SRC)/%.cpp
	g++ $(FLAGS) -c $< -o $@ $(INCLUDE)


#DEBUG PART

.PHONY: deb
deb: $(DEBUG_FOLDER) $(DEB)
	valgrind ./$(DEB) ${ARGS}

$(DEBUG_FOLDER):
	mkdir $(DEBUG_FOLDER)

$(DEB): $(DEB_C) $(DEB_CPP)
	g++ $(FLAGS) -g $^ -o $@ $(LINK)

$(DEBUG_FOLDER)/%.o: $(SRC)/%.cpp
	g++ $(FLAGS) -g -c $< -o $@ $(INCLUDE)

#CLEAN PART

clean:
	rm -rf $(BUILD_FOLDER) $(DEBUG_FOLDER)
mr_proper: clean
	rm -f $(EXEC)