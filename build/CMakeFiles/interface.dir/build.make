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

# Utility rule file for interface.

# Include the progress variables for this target.
include CMakeFiles/interface.dir/progress.make

CMakeFiles/interface: ss_tree_interface
	./ss_tree_interface

interface: CMakeFiles/interface
interface: CMakeFiles/interface.dir/build.make

.PHONY : interface

# Rule to build all files generated by this target.
CMakeFiles/interface.dir/build: interface

.PHONY : CMakeFiles/interface.dir/build

CMakeFiles/interface.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/interface.dir/cmake_clean.cmake
.PHONY : CMakeFiles/interface.dir/clean

CMakeFiles/interface.dir/depend:
	cd /mnt/c/labo-eda/code/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/labo-eda/code /mnt/c/labo-eda/code /mnt/c/labo-eda/code/build /mnt/c/labo-eda/code/build /mnt/c/labo-eda/code/build/CMakeFiles/interface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/interface.dir/depend
