/*******************************************************************************
*                                                                              *
* Author    :  Angus Johnson                                                   *
* Version   :  6.2.1                                                           *
* Date      :  31 October 2014                                                 *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2014                                         *
*                                                                              *
* License:                                                                     *
* Use, modification & distribution is subject to Boost Software License Ver 1. *
* http://www.boost.org/LICENSE_1_0.txt                                         *
*                                                                              *
* Attributions:                                                                *
* The code in this library is an extension of Bala Vatti's clipping algorithm: *
* "A generic solution to polygon clipping"                                     *
* Communications of the ACM, Vol 35, Issue 7 (July 1992) pp 56-63.             *
* http://portal.acm.org/citation.cfm?id=129906                                 *
*                                                                              *
* Computer graphics and geometric modeling: implementation and algorithms      *
* By Max K. Agoston                                                            *
* Springer; 1 edition (January 4, 2005)                                        *
* http://books.google.com/books?q=vatti+clipping+agoston                       *
*                                                                              *
* See also:                                                                    *
* "Polygon Offsetting by Computing Winding Numbers"                            *
* Paper no. DETC2005-85513 pp. 565-575                                         *
* ASME 2005 International Design Engineering Technical Conferences             *
* and Computers and Information in Engineering Conference (IDETC/CIE2005)      *
* September 24-28, 2005 , Long Beach, California, USA                          *
* http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf              *
*                                                                              *
*******************************************************************************/

#ifndef clipper_hpp
#define clipper_hpp

#define CLIPPER_VERSION "6.2.0"

//use_int32: When enabled 32bit ints are used instead of 64bit ints. This
//improve performance but coordinate values are limited to the range +/- 46340
//#define use_int32

//use_xyz: adds a Z member to IntPoint. Adds a minor cost to perfomance.
//#define use_xyz

//use_lines: Enables line clipping. Adds a very minor cost to performance.
//#define use_lines
  
//use_deprecated: Enables temporary support for the obsolete functions
//#define use_deprecated  

#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>

namespace ClipperLib {

enum ClipType { ctIntersection, ctUnion, ctDifference, ctXor };
enum PolyType { ptSubject, ptClip };
//By far the most widely used winding rules for polygon filling are
//EvenOdd & NonZero (GDI, GDI+, XLib, OpenGL, Cairo, AGG, Quartz, SVG, Gr32)
//Others rules include Positive, Negative and ABS_GTR_EQ_TWO (only in OpenGL)
//see http://glprogramming.com/red/chapter11.html
enum PolyFillType { pftEvenOdd, pftNonZero, pftPositive, pftNegative };

#ifdef use_int32
  typedef int cInt;
  static cInt const loRange = 0x7FFF;
  static cInt const hiRange = 0x7FFF;
#else
  typedef signed long long cInt;
  static cInt const loRange = 0x3FFFFFFF;
  static cInt const hiRange = 0x3FFFFFFFFFFFFFFFLL;
  typedef signed long long long64;     //used by Int128 class
  typedef unsigned long long ulong64;

#endif

struct IntPoint {
  cInt X = 0;
  cInt Y = 0;
#ifdef use_xyz
  cInt Z = 0;
  IntPoint(cInt x = 0, cInt y = 0, cInt z = 0): X(x), Y(y), Z(z) {};
#else
  explicit IntPoint(cInt x = 0, cInt y = 0): X(x), Y(y) {};
#endif

  friend inline bool operator== (const IntPoint& a, const IntPoint& b)
  {
    return a.X == b.X && a.Y == b.Y;
  }
  friend inline bool operator!= (const IntPoint& a, const IntPoint& b)
  {
    return a.X != b.X  || a.Y != b.Y; 
  }
};
//------------------------------------------------------------------------------

typedef std::vector< IntPoint > Path;
typedef std::vector< Path > Paths;

inline Path& operator <<(Path& poly, const IntPoint& p) {poly.push_back(p); return poly;}
inline Paths& operator <<(Paths& polys, const Path& p) {polys.push_back(p); return polys;}

std::ostream& operator <<(std::ostream &s, const IntPoint &p);
std::ostream& operator <<(std::ostream &s, const Path &p);
std::ostream& operator <<(std::ostream &s, const Paths &p);

struct DoublePoint
{
  double X = 0.;
  double Y = 0.;
  explicit DoublePoint(double x = 0, double y = 0) : X(x), Y(y) {}
  explicit DoublePoint(IntPoint ip) : X((double)ip.X), Y((double)ip.Y) {}
};
//------------------------------------------------------------------------------

#ifdef use_xyz
typedef void (*ZFillCallback)(IntPoint& e1bot, IntPoint& e1top, IntPoint& e2bot, IntPoint& e2top, IntPoint& pt);
#endif

enum InitOptions {ioReverseSolution = 1, ioStrictlySimple = 2, ioPreserveCollinear = 4};
enum JoinType {jtSquare, jtRound, jtMiter};
enum EndType {etClosedPolygon, etClosedLine, etOpenButt, etOpenSquare, etOpenRound};

class PolyNode;
typedef std::vector< PolyNode* > PolyNodes;

class PolyNode 
{ 
public:
    Path Contour {};
    PolyNodes Childs {};
    PolyNode* Parent = nullptr;

    PolyNode();
    virtual ~PolyNode() = default;

    PolyNode* GetNext() const;
    bool IsHole() const;
    bool IsOpen() const;
    int ChildCount() const;

private:
    unsigned Index      = 0; //node index in Parent.Childs
    bool m_IsOpen       = false;
    JoinType m_jointype = jtSquare;
    EndType m_endtype   = etClosedPolygon;

    PolyNode* GetNextSiblingUp() const;
    void AddChild(PolyNode& child);
    friend class Clipper; //to access Index
    friend class ClipperOffset; 
};

class PolyTree: public PolyNode
{ 
public:
    ~PolyTree() override { Clear(); }
    PolyNode* GetFirst() const;
    void Clear();
    int Total() const;
private:
    PolyNodes AllNodes;
    friend class Clipper; //to access AllNodes
};

bool Orientation(const Path &poly);
double Area(const Path &poly);
int PointInPolygon(const IntPoint &pt, const Path &path);

void SimplifyPolygon(const Path &in_poly, Paths &out_polys, PolyFillType fillType = pftEvenOdd);
void SimplifyPolygons(const Paths &in_polys, Paths &out_polys, PolyFillType fillType = pftEvenOdd);
void SimplifyPolygons(Paths &polys, PolyFillType fillType = pftEvenOdd);

void CleanPolygon(const Path& in_poly, Path& out_poly, double distance = 1.415);
void CleanPolygon(Path& poly, double distance = 1.415);
void CleanPolygons(const Paths& in_polys, Paths& out_polys, double distance = 1.415);
void CleanPolygons(Paths& polys, double distance = 1.415);

void MinkowskiSum(const Path& pattern, const Path& path, Paths& solution, bool pathIsClosed);
void MinkowskiSum(const Path& pattern, const Paths& paths, Paths& solution, bool pathIsClosed);
void MinkowskiDiff(const Path& poly1, const Path& poly2, Paths& solution);

void PolyTreeToPaths(const PolyTree& polytree, Paths& paths);
void ClosedPathsFromPolyTree(const PolyTree& polytree, Paths& paths);
void OpenPathsFromPolyTree(PolyTree& polytree, Paths& paths);

void ReversePath(Path& p);
void ReversePaths(Paths& p);

struct IntRect { cInt left; cInt top; cInt right; cInt bottom; };

//enums that are used internally ...
enum EdgeSide { esLeft = 1, esRight = 2};

//forward declarations (for stuff used internally) ...

enum Direction { dRightToLeft, dLeftToRight };

static int const Unassigned = -1;  //edge not currently 'owning' a solution
static int const Skip = -2;        //edge that would otherwise close a path

#define HORIZONTAL (-1.0E+40)
#define TOLERANCE (1.0e-20)
#define NEAR_ZERO(val) (((val) > -TOLERANCE) && ((val) < TOLERANCE))

struct TEdge {
  IntPoint Bot     {};
  IntPoint Curr    {};
  IntPoint Top     {};
  IntPoint Delta   {};
  double Dx        = 0.;
  PolyType PolyTyp {};
  EdgeSide Side    {};
  int WindDelta    = 0; //1 or -1 depending on winding direction
  int WindCnt      = 0;
  int WindCnt2     = 0; //winding count of the opposite polytype
  int OutIdx       = 0;
  TEdge *Next      = nullptr;
  TEdge *Prev      = nullptr;
  TEdge *NextInLML = nullptr;
  TEdge *NextInAEL = nullptr;
  TEdge *PrevInAEL = nullptr;
  TEdge *NextInSEL = nullptr;
  TEdge *PrevInSEL = nullptr;
};

struct IntersectNode {
  TEdge          *Edge1 = nullptr;
  TEdge          *Edge2 = nullptr;
  IntPoint        Pt {};
};

struct LocalMinimum {
  cInt  Y = 0;
  TEdge *LeftBound  = nullptr;
  TEdge *RightBound = nullptr;
};

struct OutPt {
  int       Idx  = -1;
  IntPoint  Pt   {};
  OutPt    *Next = nullptr;
  OutPt    *Prev = nullptr;
};

struct OutRec {
  int       Idx       = -1;
  bool      IsHole    = false;
  bool      IsOpen    = false;
  OutRec   *FirstLeft = nullptr;  //see comments in clipper.pas
  PolyNode *PolyNd    = nullptr;
  OutPt    *Pts       = nullptr;
  OutPt    *BottomPt  = nullptr;
};

struct Join {
  OutPt    *OutPt1 = nullptr;
  OutPt    *OutPt2 = nullptr;
  IntPoint  OffPt {};
};

struct LocMinSorter
{
  inline bool operator()(const LocalMinimum& locMin1, const LocalMinimum& locMin2)
  {
    return locMin2.Y < locMin1.Y;
  }
};

typedef std::vector < OutRec* > PolyOutList;
typedef std::vector < TEdge* > EdgeList;
typedef std::vector < Join* > JoinList;
typedef std::vector < IntersectNode* > IntersectList;

//------------------------------------------------------------------------------

//ClipperBase is the ancestor to the Clipper class. It should not be
//instantiated directly. This class simply abstracts the conversion of sets of
//polygon coordinates into edge objects that are stored in a LocalMinima list.
class ClipperBase
{
public:
  ClipperBase();
  virtual ~ClipperBase();
  bool AddPath(const Path &pg, PolyType PolyTyp, bool Closed);
  bool AddPaths(const Paths &ppg, PolyType PolyTyp, bool Closed);
  virtual void Clear();
  IntRect GetBounds();
  bool PreserveCollinear() { return m_PreserveCollinear; }
  void PreserveCollinear(bool value) { m_PreserveCollinear = value; }

protected:
  void DisposeLocalMinimaList();
  TEdge* AddBoundsToLML(TEdge *e, bool IsClosed);
  void PopLocalMinima();
  virtual void Reset();
  TEdge* ProcessBound(TEdge* E, bool IsClockwise);
  void DoMinimaLML(TEdge* E1, TEdge* E2, bool IsClosed);
  TEdge* DescendToMin(TEdge *&E);
  void AscendToMax(TEdge *&E, bool Appending, bool IsClosed);

  typedef std::vector<LocalMinimum> MinimaList;

  MinimaList::iterator m_CurrentLM {};
  MinimaList           m_MinimaList {};
  bool                 m_UseFullRange = false;
  EdgeList             m_edges {};
  bool                 m_PreserveCollinear = false;
  bool                 m_HasOpenPaths = false;
};
//------------------------------------------------------------------------------

class Clipper : public virtual ClipperBase
{
public:
  explicit Clipper(int initOptions = 0);
  ~Clipper() override;
  bool Execute(ClipType clipType,
    Paths &solution,
    PolyFillType subjFillType = pftEvenOdd,
    PolyFillType clipFillType = pftEvenOdd);
  bool Execute(ClipType clipType,
    PolyTree &polytree,
    PolyFillType subjFillType = pftEvenOdd,
    PolyFillType clipFillType = pftEvenOdd);
  bool ReverseSolution() {return m_ReverseOutput;};
  void ReverseSolution(bool value) {m_ReverseOutput = value;};
  bool StrictlySimple() {return m_StrictSimple;};
  void StrictlySimple(bool value) {m_StrictSimple = value;};
  //set the callback function for z value filling on intersections (otherwise Z is 0)
#ifdef use_xyz
  void ZFillFunction(ZFillCallback zFillFunc);
#endif
protected:
  void Reset() override;
  virtual bool ExecuteInternal();
private:
  typedef std::priority_queue<cInt> ScanbeamList;

  PolyOutList    m_PolyOuts      {};
  JoinList       m_Joins         {};
  JoinList       m_GhostJoins    {};
  IntersectList  m_IntersectList {};
  ClipType       m_ClipType      {};
  ScanbeamList   m_Scanbeam      {};
  TEdge         *m_ActiveEdges   = nullptr;
  TEdge         *m_SortedEdges   = nullptr;
  bool           m_ExecuteLocked = false;
  PolyFillType   m_ClipFillType  {};
  PolyFillType   m_SubjFillType  {};
  bool           m_ReverseOutput = false;
  bool           m_UsingPolyTree = false;
  bool           m_StrictSimple  = false;
#ifdef use_xyz
  ZFillCallback   m_ZFill; //custom callback 
#endif
  void SetWindingCount(TEdge& edge);
  bool IsEvenOddFillType(const TEdge& edge) const;
  bool IsEvenOddAltFillType(const TEdge& edge) const;
  void InsertScanbeam(cInt Y);
  cInt PopScanbeam();
  void InsertLocalMinimaIntoAEL(cInt botY);
  void InsertEdgeIntoAEL(TEdge *edge, TEdge* startEdge);
  void AddEdgeToSEL(TEdge *edge);
  void CopyAELToSEL();
  void DeleteFromSEL(TEdge *e);
  void DeleteFromAEL(TEdge *e);
  void UpdateEdgeIntoAEL(TEdge *&e);
  void SwapPositionsInSEL(TEdge *edge1, TEdge *edge2);
  bool IsContributing(const TEdge& edge) const;
  bool IsTopHorz(const cInt XPos);
  void SwapPositionsInAEL(TEdge *edge1, TEdge *edge2);
  void DoMaxima(TEdge *e);
  void ProcessHorizontals(bool IsTopOfScanbeam);
  void ProcessHorizontal(TEdge *horzEdge, bool isTopOfScanbeam);
  void AddLocalMaxPoly(TEdge *e1, TEdge *e2, const IntPoint &pt);
  OutPt* AddLocalMinPoly(TEdge *e1, TEdge *e2, const IntPoint &pt);
  OutRec* GetOutRec(int idx);
  void AppendPolygon(TEdge *e1, TEdge *e2);
  void IntersectEdges(TEdge *e1, TEdge *e2, IntPoint &pt);
  OutRec* CreateOutRec();
  OutPt* AddOutPt(TEdge *e, const IntPoint &pt);
  void DisposeAllOutRecs();
  void DisposeOutRec(PolyOutList::size_type index);
  bool ProcessIntersections(cInt topY);
  void BuildIntersectList(cInt topY);
  void ProcessIntersectList();
  void ProcessEdgesAtTopOfScanbeam(cInt topY);
  void BuildResult(Paths& polys);
  void BuildResult2(PolyTree& polytree);
  void SetHoleState(TEdge *e, OutRec *outrec);
  void DisposeIntersectNodes();
  bool FixupIntersectionOrder();
  void FixupOutPolygon(OutRec &outrec);
  bool IsHole(TEdge *e);
  bool FindOwnerFromSplitRecs(OutRec &outRec, OutRec *&currOrfl);
  void FixHoleLinkage(OutRec &outrec);
  void AddJoin(OutPt *op1, OutPt *op2, IntPoint offPt);
  void ClearJoins();
  void ClearGhostJoins();
  void AddGhostJoin(OutPt *op, IntPoint offPt);
  bool JoinPoints(Join *j, OutRec* outRec1, OutRec* outRec2);
  void JoinCommonEdges();
  void DoSimplePolygons();
  void FixupFirstLefts1(OutRec* OldOutRec, OutRec* NewOutRec);
  void FixupFirstLefts2(OutRec* OldOutRec, OutRec* NewOutRec);
#ifdef use_xyz
  void SetZ(IntPoint& pt, TEdge& e1, TEdge& e2);
#endif
};
//------------------------------------------------------------------------------

class ClipperOffset 
{
public:
  explicit ClipperOffset(double miterLimit = 2.0, double roundPrecision = 0.25);
  ~ClipperOffset();
  void AddPath(const Path& path, JoinType joinType, EndType endType);
  void AddPaths(const Paths& paths, JoinType joinType, EndType endType);
  void Execute(Paths& solution, double delta);
  void Execute(PolyTree& solution, double delta);
  void Clear();

  double MiterLimit = 0.;
  double ArcTolerance = 0.;

private:
  Paths m_destPolys {};
  Path m_srcPoly {};
  Path m_destPoly {};
  std::vector<DoublePoint> m_normals {};
  double m_delta = 0.;
  double m_sinA = 0.;
  double m_sin = 0.;
  double m_cos = 0.;
  double m_miterLim = 0.;
  double m_StepsPerRad = 0.;
  IntPoint m_lowest {};
  PolyNode m_polyNodes {};

  void FixOrientations();
  void DoOffset(double delta);
  void OffsetPoint(int j, int& k, JoinType jointype);
  void DoSquare(int j, int k);
  void DoMiter(int j, int k, double r);
  void DoRound(int j, int k);
};
//------------------------------------------------------------------------------

class clipperException : public std::runtime_error {
  public:
    explicit clipperException(const char *description)
      : std::runtime_error(description)
    {}
};
//------------------------------------------------------------------------------

} //ClipperLib namespace

#endif //clipper_hpp


