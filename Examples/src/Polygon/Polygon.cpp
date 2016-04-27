#include "Polygon.hpp"
#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <limits>
#include<vector>
#include<string>
#include<memory>
#include<sstream>
#include<fstream>
namespace Geometry
{
  
  double distance(Point2D const & a, Point2D const & b){
    // Not very efficient. This function should be implemented either
    // as friend or as a method, to access private members.
    Point2D tmp(a-b);
    return std::sqrt(tmp.x()*tmp.x()+
		     tmp.y()*tmp.y()
		     );
  }

  bool compare(const Point2D & A,const Point2D & B ){return (B.x()==A.x() && A.y()==B.y());}


  // ********************* BASE CLASS **********************
   AbstractPolygon::AbstractPolygon():All_Vertices(nullptr){};
    //Fatto
   AbstractPolygon::AbstractPolygon(std::vector<Point2D> & G,std::vector<unsigned int> const & Position, bool check)
     {
      All_Vertices=& G;
      Pos_Vertices=Position;
      if (check) this->checkConvexity();
     }
   //Fatto
   void AbstractPolygon::showMe(std::ostream & out)const
     {
     if (this->size()==0){
      out<< "empty polygon" << std::endl;
     }
     else {
      out<<"Vertices:"<<std::endl;
      out<<"    X     " << "   Y    "<<std::endl;
      for (unsigned int i=0; i<this->Pos_Vertices.size();i++)
       {
         out<<(*All_Vertices)[Pos_Vertices[i]].x()<<" " << (*All_Vertices)[Pos_Vertices[i]].y()<<std::endl;
       }
     }
     if(this->isconvex) std::cout<<" Polygon is convex"<<std::endl;
     else std::cout<<" Polygon is not convex"<<std::endl;
     }
   //Fatto
   void AbstractPolygon::checkConvexity()
     {
     auto mysize=this->size();
     // We consider segments and triangles as convex
     if (mysize <= 3) {
      this->isconvex=true;
      return;
     }
     //! Since we are dealing with floating points it is better to have
     //  a small number so that |a| < smallNumber means for us a==0
     double const smallNumber(1000*std::numeric_limits<double>::min());
     Point2D p;
     Point2D v;
     Point2D u;
     double res(0.0);
     double newres(0.0);
     //! C++11 sintax. decltype(expr) returns the type of the expression
     for ( decltype(mysize) i=0; i < mysize; ++i)
      {
	     p = (*All_Vertices)[Pos_Vertices[i]];
	     // ! next point
	     Point2D tmp = (*All_Vertices)[Pos_Vertices[(i+1) % mysize]];
	     v = tmp - p;
	      //! next next point
	     u = (*All_Vertices)[Pos_Vertices[(i+2) % mysize]];
	     if (i == 0) // in first loop direction is unknown, so save it in res
	           res = u.x() * v.y() - u.y() * v.x() + v.x() * p.y() - v.y() * p.x();
	     else
         {
	        newres = u.x() * v.y() - u.y() * v.x() + v.x() * p.y() - v.y() * p.x();
	        if (std::abs(res)<smallNumber){
	         // The two edges are aligned, skip test and update res
	         res=newres;
	        }
	        else if ((newres > 0 && res < 0) || (newres < 0 && res > 0) )
	         {
	           this->isconvex=false;
	           return;
	        }   
	      }
      }// end for
     this->isconvex=true;
     return;
     }
  
   Point2D & AbstractPolygon::operator [](const unsigned int P)
   {
    return (*All_Vertices)[Pos_Vertices[P]];
   }
  
  // ****   POLYGON
   //Fatto
   Polygon::Polygon(std::vector<Point2D> & G,std::vector<unsigned int> const & Position): AbstractPolygon(G,Position) {}
  
   double Polygon::area() const{
    auto siz=this->size();
    if (siz<3) return 0;
    double result(0);

    for (decltype(siz) i=0; i<siz;++i){
      // Current point I use references to avoid unnecessary construction
      // of objects
      Point2D const & p1((*All_Vertices)[Pos_Vertices[i]]);
      // Other point
      Point2D const & p2((*All_Vertices)[Pos_Vertices[(i+1) % siz]]);
      Point2D const & p0((*All_Vertices)[Pos_Vertices[i+siz-1] % siz]);
      result+=p1.x()*(p2.y()-p0.y());
    }
    return 0.5*result;
   }
  
   void Polygon::showMe(std::ostream & out)const
   {
    std::cout<<" A Generic Polygon"<<std::endl;
    AbstractPolygon::showMe(out);
   }

  // ********************* SQUARE **********************
    Square::Square(std::vector<Point2D> & G,std::vector<unsigned int> const & Position): AbstractPolygon(G,Position,false) 
     {
      this->isconvex=true;
      if(Pos_Vertices.size() != 4){
        throw std::runtime_error(" A square must be created giving four vertices");
      }
      // Check if it is a square!
      double l1=distance((*All_Vertices)[Pos_Vertices[1]],(*All_Vertices)[Pos_Vertices[0]]);
      double l2=distance((*All_Vertices)[Pos_Vertices[2]],(*All_Vertices)[Pos_Vertices[3]]);
      auto ratio=std::abs(this->area())/(l1*l2);
      if (std::abs(ratio - 1.0)>10*std::numeric_limits<double>::epsilon())
        throw std::runtime_error("Vertexes do not define a square");
      }

    //Fatto
    double Square::area() const{
    if(this->size()==0) return 0.0;
    // I want the area with sign, positive if the quad is 
    // oriented counterclockwise. So the easiest thing is to use the 
    // formula using cross product (even if it is not the most efficient choice)
    Point2D v((*All_Vertices)[Pos_Vertices[1]]-(*All_Vertices)[Pos_Vertices[0]]);
    Point2D w((*All_Vertices)[Pos_Vertices[2]]-(*All_Vertices)[Pos_Vertices[0]]);
    // area = v \times w. Positive if square counterclockwise oriented
    return v.x()*w.y()-v.y()*w.x();
    ;}
    //Fatto
    void Square::showMe(std::ostream & out) const
     {
      out<<"A Square"<<std::endl;
      AbstractPolygon::showMe(out);
     }
  
  //********************* TRIANGLE **********************
   //Fatto
   Triangle::Triangle(std::vector<Point2D> & G,std::vector<unsigned int> const & Position): AbstractPolygon(G,Position,false){
     this->isconvex=true;    
     // Check if we give 3 vertices
     // We may use assert, in this case we would disable the control
     // in the released version (-DNDEBUG). We prefer here to exit the program
     if(this->size() != 3){
      throw std::runtime_error(" A triangle must be created giving three vertices");
     }
     }
  
   double Triangle::area() const
     {
      if(this->size()==0) return 0.0;
      // I use the cross product since this is a 
      // signed area!
      Point2D v((*All_Vertices)[Pos_Vertices[1]]-(*All_Vertices)[Pos_Vertices[0]]);
      Point2D w((*All_Vertices)[Pos_Vertices[2]]-(*All_Vertices)[Pos_Vertices[0]]);
      // area = 0.5* v \times w. Positive if triangle counterclockwise oriented
      return 0.5*(v.x()*w.y()-v.y()*w.x());
     }
  
   void Triangle::showMe(std::ostream & out) const{
      out<<"A Triangle"<<std::endl;
      AbstractPolygon::showMe(out);
     }

  //**********************  GRID   ****************************//


   void Grid::Build (const std::string  & name)
     {
     std::vector<Point2D> vert;
     std::vector<std::shared_ptr<AbstractPolygon>> Pol;
     if (name.empty()) {std::cout<<"Errore\n\n";}
     else
       {
       std::ifstream f(name);

       if (f.is_open())
         {                                             //Se il file è correttamente aperto
          std::string riga,valore;
          getline(f,riga);                                       //Leggo la prima riga
          if (riga.empty())
           {                                         //Se la prima riga è vuota
             std::cout<<"Il file è vuoto"<<std::endl;
             f.close();   
           }
          else
           {
              std::stringstream ss(riga);
              unsigned int NP,NR;
    	
    	        if(getline(ss,valore,' '))
               {
                 NP=stoul(valore);
               }
              else
               {
    	         std::cout<<"Il file è vuoto"<<std::endl;
    	         f.close();                                              //Chiudo
       	       }
	            if(getline(ss,valore,' '))
               {
                 NR=stoul(valore);
               }
              else
               {
                  std::cout<<"Il file è vuoto"<<std::endl;
	                f.close();                                              //Chiudo
	             }
	  
              std::vector<Point2D> vert;
              double x,y;
              for(unsigned int i=0;i<NP;i++) //leggo i punti
	             {

	    	          if(getline(f,riga))
                   {
                     std::stringstream ss1(riga);
                     getline(ss1,valore,' '); //indica il numero dell'elemento
                     getline(ss1,valore,' '); //indica la x del punto
                     x=stod(valore);
		                 getline(ss1,valore,' '); //indica la y del punto
                     y=stod(valore);
		                 vert.push_back(Point2D (x,y));
                  }
	             }
              Points=vert;
              for(unsigned int i=0;i<NR;i++)
               {
	   	           unsigned int Tipo;
                 
                 if(getline(f,riga))
		               {
      		         std::stringstream ss2(riga);
      		          unsigned int pos;
      		          std::vector<unsigned int> vert_pol;
      		          getline(ss2,valore,' '); //indica il numero dell'elemento
      		          getline(ss2,valore,' '); //indica il tipo dell'elemento
                    Tipo=stoul(valore);
      		          while(getline(ss2,valore,' '))
      		           {
      			          pos=stoul(valore);
      			          vert_pol.push_back(pos);
      		           }
                    if(Tipo==0) 
      		           {
      			          Triangle t(Points,vert_pol); 
      			          Pol.push_back(std::make_shared<Triangle> (t));		
                     }
      		          else if (Tipo==1)
                     {
       		           Square t(Points,vert_pol);  
                     Pol.push_back(std::make_shared<Square> (t));
                     }
        		        else 
                     {   
                      Polygon t(Points,vert_pol);
      		            Pol.push_back(std::make_shared<Polygon> (t));
                     }
                    bool entra=1;
                   unsigned int cont=0;
                    for(unsigned int K=0; K<vert_pol.size();K++)
                     {

                       while (entra && cont<All_Edges.size())
                       { 
                         std::cout<<"Al solido "<<i<<" vertice "<<K<<" cont="<<cont<<" ho "<<All_Edges[cont].compared(vert_pol[K],vert_pol[(K+1)%vert_pol.size()])<<"\n\n";
                         if(All_Edges[cont].compared(vert_pol[K],vert_pol[(K+1)%vert_pol.size()]))
                          {
                            entra=0;
                          }
                         else
                          {
                            cont++;
                          }
                        
                       }
                       if (cont==All_Edges.size())
                        {
                          All_Edges.push_back(Edge(vert_pol[K],vert_pol[(K+1)%vert_pol.size()]));
                        }
                       cont=0;
                       entra=1;
                     }
                     std::cout<<"Dimensione All ="<<All_Edges.size()<<"\n\n";

                    }
               }
              f.close();
              Figure=Pol;
              for(auto K: All_Edges)
              {
                if(!(K.Edge_Intern(Figure,Points)))
                {
                  Internal_Edges.push_back(K);
                }
              }

           }    
         }
       
       }
    }

   Grid::Grid (const Grid & G)
     {
      Points=G.Points;
      Figure=G.Figure;
      All_Edges=G.All_Edges;
      Internal_Edges=G.Internal_Edges;
     }

   Grid & Grid::operator=(const Grid & G)
     {
      Points=G.Points;
      Figure=G.Figure;
      All_Edges=G.All_Edges;
      Internal_Edges=G.Internal_Edges;
      return *this;
     }

   double Grid::Full_Area()
     {
      double Area=0.;
      for(unsigned int i=0; i<Figure.size(); i++)
      {
         Area+=Figure[i]->area();
      }
      return Area;
     }

   void Grid::printG(std::ostream & out, unsigned int All) 
   {
      if(All==0)
      {
        out<<"All Edges of the Grid"<<std::endl;
        for(unsigned int i=0; i<All_Edges.size();i++)
        {
          All_Edges[i].print(out,i);
        }
      }
      else if(All==1)
      {
        out<<"All External Edges of the Grid"<<std::endl;
        for(unsigned int i=0; i<Internal_Edges.size();i++)
        {
          Internal_Edges[i].print(out,i);
        }
      }
      else
      {
        unsigned int cont=0;
        bool entra=1;
        out<<"All Internal Edges of the Grid"<<std::endl;
        for(unsigned int i=0; i<All_Edges.size();i++)
        {
          while(cont<Internal_Edges.size()  && entra)
          {
            if(All_Edges[i].compared(Internal_Edges[cont]))
              {
                entra=0;
              }
            else
              {
                cont++;
              }
          }
          if(entra==1)
           {
              All_Edges[i].print(out,i);
           }
           entra=1;
           cont=0;
        }
      }

    }

  //*********************   EDGE   ****************************


      bool Edge::Edge_Intern(const std::vector<std::shared_ptr<AbstractPolygon>> &  IPol,const std::vector<Point2D> & Points)
       {
         Point2D p1(Points[pos1]), p2(Points[pos2]);
         int cont=0;
         for(auto &i : IPol)
         {
          int dim=i->size();
          for(int j=0;j<dim;j++)
           {
             Point2D p3( (*i)[j]), p4((*i)[(j+1)%dim]);
             if ( (compare(p1,p3) && compare(p4,p2))|| (compare(p2,p3) && compare(p4,p1))){cont++;}
           }
          if(cont>1){return 1;}
         }
         return 0;
       }

     void Edge::print(std::ostream & out, unsigned int Pos) 
     {
       out<<"Edge["<<Pos<<"]= [ "<<pos1<<" "<<pos2<<" ]"<<std::endl;
     }
}
