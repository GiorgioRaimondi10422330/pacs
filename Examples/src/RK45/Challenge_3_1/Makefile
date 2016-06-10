CXXFLAGS= -std=c++11
INCLUDE= -I/home/pacs_student/Material/Examples/src/RK45/Challenge_3_1
HEADER= MyRK_1D.hpp method.hpp

.PHONY: all clean distclean exec main
.DEFAULT_GOAL = all

all: 
	@echo " "
	@echo "make main compila main_MyRK.cpp e method.cpp"
	@echo " "
	@echo "make exec compila main_MyRK.cpp e method.cpp ed esegue main_MyRK"
	@echo " "


main: main_MyRK.o method.o
	g++ main_MyRK.o method.o $(INCLUDE) -o main_MyRK

main_MyRK.o: main_MyRK.cpp $(HEADER)
	g++ main_MyRK.cpp $(CXXFLAGS) $(INCLUDE) -c -o main_MyRK.o

method.o: method.cpp $(HEADER)
	g++ method.cpp $(CXXFLAGS) $(INCLUDE) -c -o method.o
exec: main main_MyRK
	./main_MyRK

clean:
	$(RM) -f main_MyRK.o  method.o

distclean: clean
	$(RM) -f main_MyRK result.tex


