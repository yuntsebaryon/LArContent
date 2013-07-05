ifndef PROJECT_DIR
    PROJECT_DIR = $(PANDORA_DIR)/LArContent
    PROJECT_LIBRARY_DIR = $(PANDORA_DIR)/lib
else
    PROJECT_LIBRARY_DIR = $(PROJECT_DIR)/lib
endif

ifdef MONITORING
    DEFINES = -DMONITORING=1
endif

INCLUDES  = -I$(PROJECT_DIR)/include
INCLUDES += -I$(PANDORA_DIR)/PandoraSDK/include
ifdef MONITORING
    INCLUDES += -I$(PANDORA_DIR)/PandoraMonitoring/include
endif

CC = g++
CFLAGS = -c -Wall -g -w -fPIC -O2
ifdef BUILD_32BIT_COMPATIBLE
    CFLAGS += -m32
endif

SOURCES  = $(wildcard $(PROJECT_DIR)/src/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArClusterAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArClustering/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArClusterSeedAssociation/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArClusterSplitting/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArHelpers/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArMonitoring/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArObjects/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArReclustering/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArThreeDSeed/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArTwoDSeed/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArUtility/*.cc)
SOURCES += $(wildcard $(PROJECT_DIR)/src/LArVertex/*.cc)

OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(OBJECTS:.o=.d)

LIBS = -L$(PANDORA_DIR)/lib -lPandoraSDK

ifdef MONITORING
    LIBS += -lPandoraMonitoring
endif

ifdef BUILD_32BIT_COMPATIBLE
    LIBS += -m32
endif

LDFLAGS = $(LIBS) -Wl,-rpath

LIBRARY = $(PROJECT_LIBRARY_DIR)/libLArContent.so

all: $(SOURCES) $(OBJECTS)
	$(CC) $(OBJECTS) $(LIBS) -shared -o $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	$(CC) $(LDFLAGS) -fPIC $(OBJECTS) -o $@

-include $(DEPENDS)

%.o:%.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(DEFINES) -MP -MMD -MT $*.o -MT $*.d -MF $*.d -o $*.o $*.cc

install:
ifdef INCLUDE_TARGET
	rsync -r --exclude=.svn $(PROJECT_DIR)/include/ ${INCLUDE_TARGET}
endif
ifdef LIB_TARGET
	cp $(PROJECT_LIBRARY_DIR)/libLArContent.so ${LIB_TARGET}
endif

clean:
	rm -f $(OBJECTS)
	rm -f $(DEPENDS)
	rm -f $(LIBRARY)
