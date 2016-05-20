OBJS = fem1d.o gauss_seidel.o mesh.o
HEADERS = config.h

fem1d : $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@ $(LIBS)

$(OBJS) : %.o : %.cpp %.h $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -c $<

.PHONY: clean distclean

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) fem1d
