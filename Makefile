BUILD_DIR := build
INCLUDE_DIR := include
SRC_DIR := src
OPTFLAGS := -Ofast -flto -g3 -ggdb -fno-omit-frame-pointer
COMPILER_FLAGS := -std=c++14 -march=native $(OPTFLAGS) -fno-rtti -Wall -Wextra -Werror -pedantic -Wshadow -Wmissing-include-dirs -Winvalid-pch -Wformat=2
CXXFLAGS := $(COMPILER_FLAGS) -DNDEBUG -I$(INCLUDE_DIR)
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
DEP := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.d)
TARGET := $(SRC:$(SRC_DIR)/%.cpp=%)

all: $(TARGET)

$(TARGET): subdirs

subdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir $@

%: $(BUILD_DIR)/%.o
	$(CXX) $(OPTFLAGS) $< $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	$(RM) $(OBJ) $(TARGET) $(DEP)
	-@$(RM) -d $(BUILD_DIR) 2> /dev/null

.SECONDARY: $(OBJ)

.PHONY: all, clean

-include $(DEP)

