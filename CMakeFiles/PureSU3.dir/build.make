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
CMAKE_SOURCE_DIR = /home/gerrit/Documents/Projects/PureSU3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gerrit/Documents/Projects/PureSU3

# Include any dependencies generated for this target.
include CMakeFiles/PureSU3.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PureSU3.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PureSU3.dir/flags.make

CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o: src/LatticeSettings.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/LatticeSettings.cpp

CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/LatticeSettings.cpp > CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.i

CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/LatticeSettings.cpp -o CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.s

CMakeFiles/PureSU3.dir/src/Matrix.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/Matrix.cpp.o: src/Matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/PureSU3.dir/src/Matrix.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/Matrix.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/Matrix.cpp

CMakeFiles/PureSU3.dir/src/Matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/Matrix.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/Matrix.cpp > CMakeFiles/PureSU3.dir/src/Matrix.cpp.i

CMakeFiles/PureSU3.dir/src/Matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/Matrix.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/Matrix.cpp -o CMakeFiles/PureSU3.dir/src/Matrix.cpp.s

CMakeFiles/PureSU3.dir/src/Random.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/Random.cpp.o: src/Random.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/PureSU3.dir/src/Random.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/Random.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/Random.cpp

CMakeFiles/PureSU3.dir/src/Random.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/Random.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/Random.cpp > CMakeFiles/PureSU3.dir/src/Random.cpp.i

CMakeFiles/PureSU3.dir/src/Random.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/Random.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/Random.cpp -o CMakeFiles/PureSU3.dir/src/Random.cpp.s

CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o: src/SUNLattice.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/SUNLattice.cpp

CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/SUNLattice.cpp > CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.i

CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/SUNLattice.cpp -o CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.s

CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o: src/TwoLinkOperator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/TwoLinkOperator.cpp

CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/TwoLinkOperator.cpp > CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.i

CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/TwoLinkOperator.cpp -o CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.s

CMakeFiles/PureSU3.dir/src/Vectors.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/Vectors.cpp.o: src/Vectors.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/PureSU3.dir/src/Vectors.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/Vectors.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/Vectors.cpp

CMakeFiles/PureSU3.dir/src/Vectors.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/Vectors.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/Vectors.cpp > CMakeFiles/PureSU3.dir/src/Vectors.cpp.i

CMakeFiles/PureSU3.dir/src/Vectors.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/Vectors.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/Vectors.cpp -o CMakeFiles/PureSU3.dir/src/Vectors.cpp.s

CMakeFiles/PureSU3.dir/src/main.cpp.o: CMakeFiles/PureSU3.dir/flags.make
CMakeFiles/PureSU3.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/PureSU3.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PureSU3.dir/src/main.cpp.o -c /home/gerrit/Documents/Projects/PureSU3/src/main.cpp

CMakeFiles/PureSU3.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PureSU3.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/gerrit/Documents/Projects/PureSU3/src/main.cpp > CMakeFiles/PureSU3.dir/src/main.cpp.i

CMakeFiles/PureSU3.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PureSU3.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/gerrit/Documents/Projects/PureSU3/src/main.cpp -o CMakeFiles/PureSU3.dir/src/main.cpp.s

# Object files for target PureSU3
PureSU3_OBJECTS = \
"CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o" \
"CMakeFiles/PureSU3.dir/src/Matrix.cpp.o" \
"CMakeFiles/PureSU3.dir/src/Random.cpp.o" \
"CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o" \
"CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o" \
"CMakeFiles/PureSU3.dir/src/Vectors.cpp.o" \
"CMakeFiles/PureSU3.dir/src/main.cpp.o"

# External object files for target PureSU3
PureSU3_EXTERNAL_OBJECTS =

PureSU3: CMakeFiles/PureSU3.dir/src/LatticeSettings.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/Matrix.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/Random.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/SUNLattice.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/TwoLinkOperator.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/Vectors.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/src/main.cpp.o
PureSU3: CMakeFiles/PureSU3.dir/build.make
PureSU3: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
PureSU3: /usr/lib/x86_64-linux-gnu/libpthread.so
PureSU3: CMakeFiles/PureSU3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/gerrit/Documents/Projects/PureSU3/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable PureSU3"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PureSU3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PureSU3.dir/build: PureSU3

.PHONY : CMakeFiles/PureSU3.dir/build

CMakeFiles/PureSU3.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PureSU3.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PureSU3.dir/clean

CMakeFiles/PureSU3.dir/depend:
	cd /home/gerrit/Documents/Projects/PureSU3 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gerrit/Documents/Projects/PureSU3 /home/gerrit/Documents/Projects/PureSU3 /home/gerrit/Documents/Projects/PureSU3 /home/gerrit/Documents/Projects/PureSU3 /home/gerrit/Documents/Projects/PureSU3/CMakeFiles/PureSU3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/PureSU3.dir/depend
