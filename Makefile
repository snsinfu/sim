OPTFLAGS = \
  -msse4.2 \
  -O2      \
  -flto

DBGFLAGS = \
  -g       \
  -fsanitize=address

INCLUDES = \
  -I include

WARNFLAGS = \
  -Wno-c99-extensions

CXXFLAGS =          \
  -std=c++2a        \
  -pedantic         \
  -Wall             \
  -Wextra           \
  -Wconversion      \
  -Wsign-conversion \
  -Wshadow          \
  $(WARNFLAGS)      \
  $(OPTFLAGS)       \
  $(INCLUDES)

ARTIFACTS = \
  example   \
  *.dSYM

.PHONY: clean

example: example.cc include/**/*.hpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ example.cc

clean:
	rm -rf $(ARTIFACTS)
