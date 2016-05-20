OBJS = sparse.o test_sparse.o
HEADERS = sparse.h
CPPFLAGS ?= -I../sparse_matrix-0.1/ -I.
LDFLAGS ?= -L../sparse_matrix-0.1/
LIBS = 

sparse : $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $@ -lsparse_matrix $(LIBS)

$(OBJS) : %.o : %.cpp %.h $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

.PHONY: clean distclean

clean :
	$(RM) $(OBJS)

distclean : clean
	$(RM) sparse
