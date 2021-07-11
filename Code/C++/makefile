# Set project directory one level above of Makefile directory. $(CURDIR) is a GNU make variable containing the path to the current working directory
PROJDIR := $(realpath $(CURDIR))
SOURCEDIR := $(PROJDIR)/src
BUILDDIR := $(PROJDIR)/Build

# Name of the final executable
TARGET = fastConc

# Name the compiler
 CC = mpiicc -O2 -xhost -DMKL_DIRECT_CALL_SEQ -std=c++11
# CC = mpiicc -g -DMKL_DIRECT_CALL_SEQ -std=c++11

# Name the libraries
LIB_PATH = vsc-hard-mounts/leuven-apps/skylake/2018a/software/imkl/2018.1.163-iimpi-2018a/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
LIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

# Create the list of directories
DIRS = . meteorology meteorology/openField dictionary meteorology/CouetteFlow 
ifeq ($(TARGET), slowConc)
DIRS += densityKernel
else ifeq ($(TARGET), slowMISE)
DIRS += densityKernel
else ifeq ($(TARGET), fastConc)
DIRS += densityKernel
else ifeq ($(TARGET), fastMISE)
DIRS += densityKernel
endif
SOURCEDIRS = $(foreach dir, $(DIRS), $(addprefix $(SOURCEDIR)/, $(dir)))
TARGETDIRS = $(foreach dir, $(DIRS), $(addprefix $(BUILDDIR)/, $(dir)))

# Generate the GCC includes parameters by adding -I before each source folder
INCLUDES = $(foreach dir, $(SOURCEDIRS), $(addprefix -I, $(dir)))

# Add this list to VPATH, the place make will look for the source files
VPATH = $(SOURCEDIRS)

# Create a list of *.cpp sources in DIRS
SOURCES = $(foreach dir,$(SOURCEDIRS),$(wildcard $(dir)/*.cpp))

# Define objects for all sources
OBJS := $(subst $(SOURCEDIR),$(BUILDDIR),$(SOURCES:.cpp=.o))

# Define dependencies files for all objects
DEPS = $(OBJS:.o=.d)

# OS specific part
MKDIR = mkdir -p
ERRIGNORE = 2>/dev/null
SEP=/

# Remove space after separator
PSEP = $(strip $(SEP))

$(BUILDDIR)/%.o: $(SOURCEDIR)/%.cpp
	@echo Building $@
	$(CC) -MMD -c -o $@ $< $(INCLUDES) 

default: $(TARGET)

all: directories $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -L/$(LIB_PATH) -o $(TARGET) $^ $(LIB)

# Include dependencies
 -include $(DEPS)

directories: 
	$(MKDIR) $(subst /,$(PSEP),$(TARGETDIRS)) $(ERRIGNORE)

clean:
	rm -r Build; rm $(TARGET);
