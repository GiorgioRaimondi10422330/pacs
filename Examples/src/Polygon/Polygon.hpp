#ifndef HH_POLYGON_HH
#define HH_POLYGON_HH
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include<vector>
#include<string>
#include<sstream>
#include<memory>
#include<fstream>
/*!  @file Polygon.hpp 
  @brief This is an example of class hierarchy that
  represents poligonal object.

  @detail This module is wrapped in the namespace `Geometry`. It
  represent a first example of class hyerarchy with public
  inheritance. The class `AbstractPolygon` defines the general public
  interface of all other classes representing polygonal objects.
  
  We make use of some syntax of the C++11 new standard so you need to
  compile with *-std=c++11*
   */
namespace Geometry
{
  //! A class that holds 2D points
    /*! It also represents a vector in R2
   */
  class Point2D
   {
   public:
    //! Constructor giving coordinates.
    Point2D(double xx=0.0, double yy=0.0):coor{{xx,yy}}{}
    //! Copy constructor
    Point2D(const Point2D&)=default;
    //! Returns coordinates in a array<double>.
    std::array<double,2> get() const { return coor;}
    //! Sets point coordinates
    void set(double const &xx, double const &yy)
     {
      coor={{xx,yy}};
     }
    //! x coordinate
    double x() const {return coor[0];}
    //! y coordinate
    double y() const {return coor[1];}

    //! Subtraction is implemented as an external friend function
    friend Point2D operator - (Point2D const & a, Point2D const & b);
    //! Addition is implemented as an external friend function
    friend Point2D operator + (Point2D const & a, Point2D const & b);
   private:
    std::array<double,2> coor;
   };

  //! An alias
  using R2Vector=Point2D;

  //! subtraction operator.
  inline Point2D operator - (Point2D const & a, Point2D const & b){
    return Point2D(a.coor[0]-b.coor[0],
                   a.coor[1]-b.coor[1]);
   }

   //! Addition operator.
   /*!  
    It is defined in the header file because I want it to be
    inline.
  */
  inline Point2D operator + (Point2D const & a, Point2D const & b){
    return Point2D(a.coor[0]+b.coor[0],
                   a.coor[1]+b.coor[1]);
  }

  bool compare(const Point2D & A,const Point2D & B );
  //! Distance between points
  double distance(Point2D const & a, Point2D const & b);
 
  //! Polygon vertices are just vectors of points.
  using Vertices=std::vector<Point2D>;

  
  //Fatto
  class AbstractPolygon
   {
   public:
    //Fatto
    AbstractPolygon(std::vector<Point2D> & G, std::vector<unsigned int> const & Position ,bool check=true);
    //Fatto
    AbstractPolygon();
    //Fatto
    AbstractPolygon & operator=(AbstractPolygon const&)=default;
    //Fatto
    AbstractPolygon(AbstractPolygon const &)=default;
    //Fatto
    AbstractPolygon(AbstractPolygon&&)=default;
    //! Fatto
    AbstractPolygon & operator=(AbstractPolygon&&)=default;
    // Fatto... NB non chiamo distuctor per il puntatore perchè appartiene a Grid
    virtual ~AbstractPolygon(){};
    //! Returns the number of vertices.
    /*!  We return Vertices::size_type and not just int because
      size_type is guaranteed to be the correct type for indexes in
      stl vectors. Its actual type may be implementation dependent.
      In this case, however, int would have been fine (size_type is
      guaranteed to be an integral type, more precisely
      a type convertible to unsigned int).
    */
    //Fatto
    Vertices::size_type size() const {return Pos_Vertices.size();}
    //Fatto
    bool isConvex() const {return isconvex;}

    Point2D & operator [](const unsigned int P);
    //Fatto
    virtual void showMe(std::ostream & out=std::cout) const;
    //Fatto
    virtual double area() const=0;
   protected:
    std::vector<Point2D>* All_Vertices;
    std::vector<unsigned int> Pos_Vertices;

    bool isconvex;
    //! Test convexity of the polygon
    void checkConvexity();
   };
  //Fatto
  class Polygon: public AbstractPolygon
   {
   public:
    //
    //Fatto;
    Polygon(std::vector<Point2D> & G,std::vector<unsigned int> const & Position);
    //Fatto
    virtual ~Polygon(){};
    //Fatto
    virtual double area() const;
    //Fatto
    virtual void showMe(std::ostream & out=std::cout) const;
   };

  class Square final: public AbstractPolygon
   {
   public:
    //Fatto
    Square(std::vector<Point2D> & G,std::vector<unsigned int> const & Position);
    //Fatto
    Square(Square const &)=default;
    //Fatto
    Square(Square&&)=default;
    //Fatto
    Square & operator=(const Square &)=default;
    //Fatto
    Square & operator=(Square &&)=default;
    //! Specialised version for squares
    double area() const;
    //! Specialised version for squares.
    void showMe(std::ostream & out=std::cout) const;
   }; 
  //! A triangle
  class Triangle final: public AbstractPolygon
   {
   public:
    //Fatto
    Triangle( std::vector<Point2D> & G,std::vector<unsigned int> const & Position);
    //Fatto
    Triangle(Triangle const &)=default;
    //Fatto
    Triangle(Triangle&&)=default;
    //Fatto
    Triangle & operator=(const Triangle &)=default;
    //Fatto
    Triangle & operator=(Triangle &&)=default;
    //Fatto
    virtual double area() const;
    //Fatto
    virtual void showMe(std::ostream & out=std::cout) const;
   };

  class Edge
   {
    unsigned int pos1, pos2;

    public:

    Edge():pos1(0),pos2(0){};

    Edge(const Edge & E):pos1(E.pos1),pos2(E.pos2){};

    Edge(const unsigned int p1, const unsigned int p2):pos1(p1),pos2(p2){};

    Edge operator =(const Edge E){ pos1=E.pos1; pos2=E.pos2; return *this;};

    bool Edge_Intern(const std::vector<std::shared_ptr<AbstractPolygon>> &  IPol,const std::vector<Point2D> & Points);

    bool compared(unsigned int p1, unsigned int p2 )    {      return ( (pos1==p1 && pos2==p2 ) || (pos1==p2 && pos2==p1 ) );    }

    bool compared(const Edge & E){return ( (pos1==E.pos1 && pos2==E.pos2 ) || (pos1==E.pos2 && pos2==E.pos1 ) );}

    void print(std::ostream & out, unsigned int Pos=0) ;
   };
  
  class Grid
   {
     std::vector<Point2D> Points;
     std::vector<std::shared_ptr<AbstractPolygon> > Figure;
     std::vector<Edge> All_Edges, Internal_Edges;
     friend class Triangle;
     friend class Square;
     friend class Polygon;
     friend class Edge;
     friend class AbstractPolygon;
     
     public:
     
     Grid ()=default;
     void Build(const std::string  & name);
     Grid(const Grid & G);
     Grid & operator=(const Grid & G);
     double Full_Area();
     void printG(std::ostream & out, unsigned int All=0) ;
   };
  



};

#endif
