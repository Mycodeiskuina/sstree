# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/labo-eda/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/labo-eda/code/build

# Include any dependencies generated for this target.
include CMakeFiles/ss_tree_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ss_tree_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ss_tree_test.dir/flags.make

CMakeFiles/ss_tree_test.dir/main.cpp.o: CMakeFiles/ss_tree_test.dir/flags.make
CMakeFiles/ss_tree_test.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/labo-eda/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ss_tree_test.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ss_tree_test.dir/main.cpp.o -c /mnt/c/labo-eda/code/main.cpp

CMakeFiles/ss_tree_test.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ss_tree_test.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/labo-eda/code/main.cpp > CMakeFiles/ss_tree_test.dir/main.cpp.i

CMakeFiles/ss_tree_test.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ss_tree_test.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/labo-eda/code/main.cpp -o CMakeFiles/ss_tree_test.dir/main.cpp.s

CMakeFiles/ss_tree_test.dir/SStree.cpp.o: CMakeFiles/ss_tree_test.dir/flags.make
CMakeFiles/ss_tree_test.dir/SStree.cpp.o: ../SStree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/labo-eda/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ss_tree_test.dir/SStree.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ss_tree_test.dir/SStree.cpp.o -c /mnt/c/labo-eda/code/SStree.cpp

CMakeFiles/ss_tree_test.dir/SStree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ss_tree_test.dir/SStree.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/labo-eda/code/SStree.cpp > CMakeFiles/ss_tree_test.dir/SStree.cpp.i

CMakeFiles/ss_tree_test.dir/SStree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ss_tree_test.dir/SStree.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/labo-eda/code/SStree.cpp -o CMakeFiles/ss_tree_test.dir/SStree.cpp.s

# Object files for target ss_tree_test
ss_tree_test_OBJECTS = \
"CMakeFiles/ss_tree_test.dir/main.cpp.o" \
"CMakeFiles/ss_tree_test.dir/SStree.cpp.o"

# External object files for target ss_tree_test
ss_tree_test_EXTERNAL_OBJECTS =

ss_tree_test: CMakeFiles/ss_tree_test.dir/main.cpp.o
ss_tree_test: CMakeFiles/ss_tree_test.dir/SStree.cpp.o
ss_tree_test: CMakeFiles/ss_tree_test.dir/build.make
ss_tree_test: CMakeFiles/ss_tree_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/labo-eda/code/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ss_tree_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ss_tree_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ss_tree_test.dir/build: ss_tree_test

.PHONY : CMakeFiles/ss_tree_test.dir/build

CMakeFiles/ss_tree_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ss_tree_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ss_tree_test.dir/clean

CMakeFiles/ss_tree_test.dir/depend:
	cd /mnt/c/labo-eda/code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/labo-eda/code /mnt/c/labo-eda/code /mnt/c/labo-eda/code/build /mnt/c/labo-eda/code/build /mnt/c/labo-eda/code/build/CMakeFiles/ss_tree_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ss_tree_test.dir/depend

