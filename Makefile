.PHONY: clean

SRC_DIR = src
INC_DIR = include
LIB_DIR = lib

SOURCES = $(SRC_DIR)/main.cpp $(SRC_DIR)/geometry.cpp
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = program.exe

CXXFLAGS = -I$(INC_DIR) 
LDFLAGS = -L$(LIB_DIR) -lccd

$(TARGET): $(OBJECTS)
	g++ $^ -o $@ $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(SRC_DIR)/*.o $(TARGET)