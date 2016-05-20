CXXFLAGS=-std=c++11
CPPFLAGS= -I$(mkEigenInc) -I/home/pacs_student/usr/include -I$(PWD)
LIBS= -lmuparser -ldl -ladaptive_quadrature -lmidpoint -ltrapezoidal
LIBTEST=-ldl -ladaptive_quadrature -lmidpoint -ltrapezoidal
LDFLAGS= -L/home/pacs_student/usr/lib -L$(PWD)

.PHONY : clean distclean all

all : adaptive_quadrature.so midpoint.so trapezoidal.so fem1d test_adaptive_quadrature 

midpoint.so : midpoint.cpp midpoint.hpp
	$(CXX) $(CPPFLAGS) -fPIC -shared $(CXXFLAGS) midpoint.cpp -o libmidpoint.so

adaptive_quadrature.so : adaptive_quadrature.cpp adaptive_quadrature.hpp
	$(CXX) $(CPPFLAGS) -fPIC -shared $(CXXFLAGS) adaptive_quadrature.cpp -o libadaptive_quadrature.so

trapezoidal.so : trapezoidal.cpp trapezoidal.hpp
	$(CXX) $(CPPFLAGS) -fPIC -shared $(CXXFLAGS) trapezoidal.cpp -o libtrapezoidal.so

test_adaptive_quadrature.o : test_adaptive_quadrature.cpp 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c test_adaptive_quadrature.cpp

test_adaptive_quadrature : test_adaptive_quadrature.o adaptive_quadrature.so
	$(CXX) $(LDFLAGS) $(CPPFLAGS) -o test_adaptive_quadrature test_adaptive_quadrature.o $(LIBTEST)

fem1d : fem1d.o mesh.o 
	$(CXX) $(LDFLAGS) fem1d.o mesh.o -o $@ $(LIBS)

fem1d.o mesh.o : %.o : %.cpp %.h $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(LIBS) $<

clean :
	$(RM) adaptive_quadrature.o test_adaptive_quadrature.o fem1d.o mesh.o midpoint.o trapezoidal.o 

distclean : clean
	$(RM) test_adaptive_quadrature adaptive_quadrature.so midpoint.so trapezoidal.so fem1d
