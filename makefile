
SOURCES = _main.cpp \
    bdmpx.c \
	csv.c \
	dotproduct.cpp \
    rsqrt.cpp \
	sndmix.cpp \
	sndpaint.cpp \
	diffusecolor.cpp \
	lerpmeshvertexes.cpp \
	projectdlighttexture.cpp

#DEFINES = -mfloat-abi=softfp -mfpu=neon -fno-default-inline -Winline -fpermissive -funsafe-math-optimizations
DEFINES = -fno-default-inline -Winline -fpermissive -O1

INCLUDES = -I$(STAGING_INC) \
            -I. \
            -I./include

LIB_PATH = -L.

LIBS = -lrt

TARGET = testdata/test

SUBUNITS = 

INSTALL = install

OBJECTS_I = $(SOURCES:%.c=obj/arm/%.o)
OBJECTS = $(OBJECTS_I:%.cpp=obj/arm/%.o)


all: $(OBJECTS)
	@list='$(SUBUNITS)'; \
	for subunit in $$list; do \
		echo "Building $$subunit"; \
		(cd $$subunit && $(MAKE) all) \
	done;
	$(CXX) $(CXXFLAGS) $(DEFINES) -o $(TARGET) $(OBJECTS) $(LDFLAGS) $(LIB_PATH) $(LIBS)

# make sure sources and dirs exist before compiling objects
$(OBJECTS): $(SOURCES) makedirs

obj/arm/%.o: %.c
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -save-temps -c $< -o $@

obj/arm/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDES) -save-temps -c $< -o $@

makedirs:
	@if [ ! -d "obj" ]; then mkdir -p obj; fi
	@if [ ! -d "obj/arm" ]; then mkdir -p obj/arm; fi

clean:
	@list='$(SUBUNITS)'; \
	for subunit in $$list; do \
		echo "Cleaning $$subunit"; \
		(cd $$subunit && $(MAKE) clean) \
	done;
	rm -rf $(TARGET) $(OBJECTS) *.ii *.s
	rm -rf ./obj

install:
	$(INSTALL) --mode=755 -dv $(DESTDIR)$(INSTALL_BINDIR)
	$(INSTALL) --mode=755 $(TARGET) $(DESTDIR)$(INSTALL_BINDIR)

