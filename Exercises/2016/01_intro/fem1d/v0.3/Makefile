OBJS = fem1d.o
HEADERS = fem1d.h mesh.h gauss_seidel.h config.h

fem1d : fem1d.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -c $<
	$(CXX) $(LDFLAGS) -std=c++11 $@.o -o $@ $(LIBS)

.PHONY: clean distclean

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) fem1d
