# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/amine/workspace/yasimSBML

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/amine/workspace/yasimSBML

# Include any dependencies generated for this target.
include CMakeFiles/yasimSBML.exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/yasimSBML.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/yasimSBML.exe.dir/flags.make

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o: CMakeFiles/yasimSBML.exe.dir/flags.make
CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o: src/yasimSBML.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/amine/workspace/yasimSBML/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o   -c /home/amine/workspace/yasimSBML/src/yasimSBML.c

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/amine/workspace/yasimSBML/src/yasimSBML.c > CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.i

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/amine/workspace/yasimSBML/src/yasimSBML.c -o CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.s

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.requires:
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.requires

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.provides: CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.requires
	$(MAKE) -f CMakeFiles/yasimSBML.exe.dir/build.make CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.provides.build
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.provides

CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.provides.build: CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.provides.build

CMakeFiles/yasimSBML.exe.dir/src/especes.c.o: CMakeFiles/yasimSBML.exe.dir/flags.make
CMakeFiles/yasimSBML.exe.dir/src/especes.c.o: src/especes.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/amine/workspace/yasimSBML/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/yasimSBML.exe.dir/src/especes.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/yasimSBML.exe.dir/src/especes.c.o   -c /home/amine/workspace/yasimSBML/src/especes.c

CMakeFiles/yasimSBML.exe.dir/src/especes.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/yasimSBML.exe.dir/src/especes.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/amine/workspace/yasimSBML/src/especes.c > CMakeFiles/yasimSBML.exe.dir/src/especes.c.i

CMakeFiles/yasimSBML.exe.dir/src/especes.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/yasimSBML.exe.dir/src/especes.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/amine/workspace/yasimSBML/src/especes.c -o CMakeFiles/yasimSBML.exe.dir/src/especes.c.s

CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.requires:
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.requires

CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.provides: CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.requires
	$(MAKE) -f CMakeFiles/yasimSBML.exe.dir/build.make CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.provides.build
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.provides

CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.provides.build: CMakeFiles/yasimSBML.exe.dir/src/especes.c.o
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.provides.build

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o: CMakeFiles/yasimSBML.exe.dir/flags.make
CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o: src/simulation.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/amine/workspace/yasimSBML/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o   -c /home/amine/workspace/yasimSBML/src/simulation.c

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/yasimSBML.exe.dir/src/simulation.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/amine/workspace/yasimSBML/src/simulation.c > CMakeFiles/yasimSBML.exe.dir/src/simulation.c.i

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/yasimSBML.exe.dir/src/simulation.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/amine/workspace/yasimSBML/src/simulation.c -o CMakeFiles/yasimSBML.exe.dir/src/simulation.c.s

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.requires:
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.requires

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.provides: CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.requires
	$(MAKE) -f CMakeFiles/yasimSBML.exe.dir/build.make CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.provides.build
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.provides

CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.provides.build: CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o
.PHONY : CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.provides.build

# Object files for target yasimSBML.exe
yasimSBML_exe_OBJECTS = \
"CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o" \
"CMakeFiles/yasimSBML.exe.dir/src/especes.c.o" \
"CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o"

# External object files for target yasimSBML.exe
yasimSBML_exe_EXTERNAL_OBJECTS =

bin/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o
bin/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/especes.c.o
bin/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o
bin/yasimSBML.exe: /usr/lib64/libxml2.so
bin/yasimSBML.exe: /usr/local/lib/libsbml.so
bin/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/build.make
bin/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable bin/yasimSBML.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/yasimSBML.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/yasimSBML.exe.dir/build: bin/yasimSBML.exe
.PHONY : CMakeFiles/yasimSBML.exe.dir/build

# Object files for target yasimSBML.exe
yasimSBML_exe_OBJECTS = \
"CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o" \
"CMakeFiles/yasimSBML.exe.dir/src/especes.c.o" \
"CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o"

# External object files for target yasimSBML.exe
yasimSBML_exe_EXTERNAL_OBJECTS =

CMakeFiles/CMakeRelink.dir/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/especes.c.o
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: /usr/lib64/libxml2.so
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: /usr/local/lib/libsbml.so
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/build.make
CMakeFiles/CMakeRelink.dir/yasimSBML.exe: CMakeFiles/yasimSBML.exe.dir/relink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable CMakeFiles/CMakeRelink.dir/yasimSBML.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/yasimSBML.exe.dir/relink.txt --verbose=$(VERBOSE)

# Rule to relink during preinstall.
CMakeFiles/yasimSBML.exe.dir/preinstall: CMakeFiles/CMakeRelink.dir/yasimSBML.exe
.PHONY : CMakeFiles/yasimSBML.exe.dir/preinstall

CMakeFiles/yasimSBML.exe.dir/requires: CMakeFiles/yasimSBML.exe.dir/src/yasimSBML.c.o.requires
CMakeFiles/yasimSBML.exe.dir/requires: CMakeFiles/yasimSBML.exe.dir/src/especes.c.o.requires
CMakeFiles/yasimSBML.exe.dir/requires: CMakeFiles/yasimSBML.exe.dir/src/simulation.c.o.requires
.PHONY : CMakeFiles/yasimSBML.exe.dir/requires

CMakeFiles/yasimSBML.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/yasimSBML.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/yasimSBML.exe.dir/clean

CMakeFiles/yasimSBML.exe.dir/depend:
	cd /home/amine/workspace/yasimSBML && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/amine/workspace/yasimSBML /home/amine/workspace/yasimSBML /home/amine/workspace/yasimSBML /home/amine/workspace/yasimSBML /home/amine/workspace/yasimSBML/CMakeFiles/yasimSBML.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/yasimSBML.exe.dir/depend

