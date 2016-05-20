OBJS = fem1d.o mesh.o
HEADERS = 
CPPFLAGS = -I/home/pacs_student/Material/Exercises/2016/02_getpot/matrix-0.1 -I. -I$(mkOpenblasLib)
LDFLAGS = -L/home/pacs_student/Material/Exercises/2016/02_getpot/matrix-0.1 
LDFLAGS+= -L$(mkOpenblasLib)
LIBS += -lopenblas
CXXFLAGS+= -std=c++11
fem1d : $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@ -lmatrix $(LIBS)

$(OBJS) : %.o : %.cpp %.h $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

.PHONY: clean distclean

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) fem1d
