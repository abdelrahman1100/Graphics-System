#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

using namespace std;

#include <iostream>
#include <tchar.h>
#include <windows.h>
#include<cmath>
#include<list>
#include<stack>
#include<vector>
#include <fstream>

#define ID_CLEAR 1001
#define IDR_MENU 1000
#define MAXENTRIES 600
#define MAXINT 1e5
const int OO = 1e9;
const int mxEntries = 600;
double const PI = 3.14159265358979323846;

//-------------------------------------------------------------------------------------------------------------------------
// MARIAM
int Round(double val) {
    return (int) val + 0.5;
}

void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color) {
    SetPixel(hdc, xc + a, yc + b, color);
    SetPixel(hdc, xc - a, yc + b, color);
    SetPixel(hdc, xc - a, yc - b, color);
    SetPixel(hdc, xc + a, yc - b, color);
    SetPixel(hdc, xc + b, yc + a, color);
    SetPixel(hdc, xc - b, yc + a, color);
    SetPixel(hdc, xc - b, yc - a, color);
    SetPixel(hdc, xc + b, yc - a, color);
}

void CircleDirect(HDC hdc, int xc, int yc, int R, COLORREF color, int border) {
    int x = 0, y = R, cnt = 0;
    int R2 = R * R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                x += 3;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                x += 5;
                continue;
            }
        }
        x++;
        y = round(sqrt((double) (R2 - x * x)));
        Draw8Points(hdc, xc, yc, x, y, color);
        cnt++;
    }
}

void CirclePolar(HDC hdc, int xc, int yc, int R, COLORREF color, int border) {
    int x = R, y = 0, cnt = 0;
    double theta = 0, dtheta = 1.0 / R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x > y) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                theta += 3 * dtheta;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                theta += 5 * dtheta;
                continue;
            }
        }
        theta += dtheta;
        x = round(R * cos(theta));
        y = round(R * sin(theta));
        Draw8Points(hdc, xc, yc, x, y, color);
        cnt++;
    }
}

void CircleIterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color, int border) {
    double x = R, y = 0;
    int cnt = 0;
    double dtheta = 1.0 / R;
    double cdtheta = cos(dtheta), sdtheta = sin(dtheta);
    Draw8Points(hdc, xc, yc, R, 0, color);
    while (x > y) {
        if (border == 2) // dotted border
        {
            if (cnt % 3 == 0) {
                cnt++;
                x = x * cdtheta - y * sdtheta;
                y = x * sdtheta + y * cdtheta;
                continue;
            }
        } else if (border == 3) // dashed border
        {
            if (cnt % 10 == 0) {
                cnt++;
                x = x * cdtheta - y * sdtheta;
                y = x * sdtheta + y * cdtheta;
                continue;
            }
        }
        double x1 = x * cdtheta - y * sdtheta;
        y = x * sdtheta + y * cdtheta;
        x = x1;
        Draw8Points(hdc, xc, yc, round(x), round(y), color);
        cnt++;
    }
}

void CircleMidPoint(HDC hdc, int xc, int yc, int R, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0)
            d += 2 * x + 2;
        else {

            d += 2 * (x - y) + 5;
            y--;
        }
        x++;
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

void CircleModifyMid(HDC hdc, int xc, int yc, int R, COLORREF color) {
    int x = 0, y = R;
    int d = 1 - R;
    int c1 = 3, c2 = 5 - 2 * R;
    Draw8Points(hdc, xc, yc, x, y, color);
    while (x < y) {
        if (d < 0) {
            d += c1;
            c2 += 2;
        } else {

            d += c2;
            c2 += 4;
            y--;
        }
        c1 += 2;
        x++;
        Draw8Points(hdc, xc, yc, x, y, color);
    }
}

//void CircleMidPoint(HDC hdc, int xc, int yc, int R, COLORREF color, int border)
//{
//    int x = 0, y = R, cnt = 0;
//    int d = 1 - R;
//    Draw8Points(hdc, xc, yc, x, y, color);
//    while (x <= y)
//    {
//        if (border == 2) // dotted border
//        {
//            if (cnt % 3 == 0)
//            {
//                cnt++;
//                x++;
//                d += 4 * x + 6;
//                continue;
//            }
//        }
//        else if (border == 3) // dashed border
//        {
//            if (cnt % 10 == 0)
//            {
//                cnt++;
//                x++;
//                d += 4 * x +6;
//                continue;
//            }
//        }
//        if (d < 0)
//        {
//            d += 4 * x + 6;
//        }
//        else
//        {
//            d += 4 * (x - y) + 10;
//            y--;
//        }
//        x++;
//        Draw8Points(hdc, xc, yc, x, y, color);
//        cnt++;
//    }
//}


//void CircleModifyMid(HDC hdc, int xc, int yc, int R, COLORREF color, int border)
//{
//    int x = 0, y = R;
//    int d = 1 - R;
//    int cnt = 0;
//    Draw8Points(hdc, xc, yc, x, y, color);
//
//    while (x <= y)
//    {
//        if (border == 2) // dotted border
//        {
//            if (cnt % 3 == 0)
//            {
//                cnt++;
//                x++;
//                d += 2 * x + 1;
//                continue;
//            }
//        }
//        else if (border == 3) // dashed border
//        {
//            if (cnt % 10 == 0)
//            {
//                cnt++;
//                x++;
//                d += 2 * x + 1;
//                continue;
//            }
//        }
//
//        if (d < 0)
//        {
//            d += 2 * x + 3;
//        }
//        else
//        {
//            d += 2 * (x - y) + 5;
//            y--;
//        }
//        x++;
//
//        Draw8Points(hdc, xc, yc, x, y, color);
//        cnt++;
//    }
//}

/////////////////////////////////////////////////////////////convex fill///////////////////////////////////////////////////////////////////////////

struct Entry {
    int xmin, xmax;
};

void InitEntries(Entry table[]) {
    for (int i = 0; i < MAXENTRIES; i++) {

        table[i].xmin = MAXINT;
        table[i].xmax = -MAXINT;

    }
}

void ScanEdge(POINT v1, POINT v2, Entry table[]) {
    if (v1.y == v2.y)return;
    if (v1.y > v2.y)swap(v1, v2);
    double minv = (double) (v2.x - v1.x) / (v2.y - v1.y);
    double x = v1.x;
    int y = v1.y;
    while (y < v2.y) {
        if (x < table[y].xmin)table[y].xmin = (int) ceil(x);
        if (x > table[y].xmax)table[y].xmax = (int) floor(x);
        y++;
        x += minv;
    }
}

void DrawSanLines(HDC hdc, Entry table[], COLORREF color) {
    for (int y = 0; y < MAXENTRIES; y++)
        if (table[y].xmin < table[y].xmax)
            for (int x = table[y].xmin; x <= table[y].xmax; x++)
                SetPixel(hdc, x, y, color);

}

void ConvexFill(HDC hdc, POINT p[], int n, COLORREF color) {
    Entry *table = new Entry[MAXENTRIES];
    InitEntries(table);
    POINT v1 = p[n - 1];
    for (int i = 0; i < n; i++) {
        POINT v2 = p[i];
        ScanEdge(v1, v2, table);
        v1 = p[i];
    }
    DrawSanLines(hdc, table, color);
    delete table;
}

//////////////////////////////////////////////non convex fill/////////////////////////////////////////////////////////////////////////////////////////////
struct EdgeRec {
    double x;
    double minv;
    int ymax;

    bool operator<(EdgeRec r) {
        return x < r.x;
    }
};

typedef list<EdgeRec> EdgeList;

EdgeRec InitEdgeRec(POINT &v1, POINT &v2) {
    if (v1.y > v2.y)swap(v1, v2);
    EdgeRec rec;
    rec.x = v1.x;
    rec.ymax = v2.y;
    rec.minv = (double) (v2.x - v1.x) / (v2.y - v1.y);
    return rec;
}

void InitEdgeTable(POINT *polygon, int n, EdgeList table[]) {
    POINT v1 = polygon[n - 1];
    for (int i = 0; i < n; i++) {
        POINT v2 = polygon[i];
        if (v1.y == v2.y) {
            v1 = v2;
            continue;
        }
        EdgeRec rec = InitEdgeRec(v1, v2);
        table[v1.y].push_back(rec);
        v1 = polygon[i];
    }
}

void GeneralPolygonFill(HDC hdc, POINT *polygon, int n, COLORREF c) {
    EdgeList *table = new EdgeList[MAXENTRIES];
    InitEdgeTable(polygon, n, table);
    int y = 0;
    while (y < MAXENTRIES && table[y].size() == 0)y++;
    if (y == MAXENTRIES)return;
    EdgeList ActiveList = table[y];
    while (ActiveList.size() > 0) {
        ActiveList.sort();
        for (EdgeList::iterator it = ActiveList.begin(); it != ActiveList.end(); it++) {
            int x1 = (int) ceil(it->x);
            it++;
            int x2 = (int) floor(it->x);
            for (int x = x1; x <= x2; x++)SetPixel(hdc, x, y, c);
        }
        y++;
        EdgeList::iterator it = ActiveList.begin();
        while (it != ActiveList.end())
            if (y == it->ymax) it = ActiveList.erase(it); else it++;
        for (EdgeList::iterator it = ActiveList.begin(); it != ActiveList.end(); it++)
            it->x += it->minv;
        ActiveList.insert(ActiveList.end(), table[y].begin(), table[y].end());
    }
    delete[] table;
}

///////////////////////////////////////////////////////////////////Ellipse////////////////////////////////////////////////////////////////////////////////
void Draw4Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c) {

    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc - y, c);
}
void DrawDirectEllipse(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {

    int x = 0;
    double y = B;
    Draw4Points(hdc, xc, yc, 0, B, c);

    while (x * B * B < y * A * A) {
        x++;
        y = B * sqrt(1 - (double) x * x / (A * A));
        Draw4Points(hdc, xc, yc, x, Round(y), c);
    }
    int x1 = A;
    double y1 = 0;
    Draw4Points(hdc, xc, yc, A, 0, c);

    while (x1 * B * B > y1 * A * A) {
        y1++;
        x1 = A * sqrt(1 - (double) y1 * y1 / (B * B));
        Draw4Points(hdc, xc, yc, Round(x1), y1, c);
    }

}
void DrawPolarEllipse(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {

    double theta = 1.0 / max(A, B), x = 0, y = A;
    double st = sin(theta);
    double ct = cos(theta);

    while (x < y) {

        double x1 = x * ct - (double) A / B * y * st;
        y = (double) B / A * x * st + y * ct;
        x = x1;
        Draw4Points(hdc, xc, yc, Round(x), Round(y), c);
    }
    while (x > y) {
        double x1 = x * ct - (double) A / B * y * st;
        y = (double) B / A * x * st + y * ct;
        x = x1;
        Draw4Points(hdc, xc, yc, Round(x), Round(y), c);
    }
}
void EllipseBresenham(HDC hdc, int xc, int yc, int A, int B, COLORREF c){
    int x = 0, y = B;
    int d1 = B*B - A*A*B + A*A/1;
    Draw4Points(hdc, xc, yc, x, y, c);
    while (A*A*(y-0.5) > B*B*(x+1)){
        if (d1 < 0){
            d1 += B*B*(2*x+3);
            x++;
        }
        else{
            d1 += B*B*(2*x+3) + A*A*(-2*y+2);
            x++;
            y--;
        }
        Draw4Points(hdc, xc, yc, x, y, c);
    }
    int d2 = B*B*(x+0.5)*(x+0.5) + A*A*(y-1)*(y-1) - A*A*B*B;
    while (y > 0){
        if (d2 < 0){
            d2 += B*B*(2*x+2) + A*A*(-2*y+3);
            x++;
            y--;
        }
        else{
            d2 += A*A*(-2*y+3);
            y--;
        }
        Draw4Points(hdc, xc, yc, x, y, c);
    }
}

//void Ellipsederict(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
//    int x = A;
//    int y = 0;
//    double dtheta = 1.0 / max(A, B);
//    double theta = 0;
//    draw4Point(hdc, xc, yc, x, y, c);
//    while (x>0) {
//        x = A * cos(theta);
//        y = B * sin(theta);
//        theta += dtheta;
//        draw4Point(hdc, xc, yc, x, y, c);
//    }
//}
//void Ellipsederict2(HDC hdc, int xc, int yc, int A, int B, COLORREF c) {
//double dtheta = 1.0 / max(A, B);
//int x=A,y=0;
//int x1=x-xc,y1=y-yc;
//int x2=x1/A;
//int y2=y1/B;
//for(double th = 0;th <= 2*PI;th += dtheta){
//int x3=x2*cos(dtheta)-y2* sin(dtheta);
//int y3=y2*cos(dtheta)+x2* sin(dtheta);
//SetPixel(hdc,x3,y3,c);
//}
//}
//------------------------------------------------------------------------------------------------------------------------------------------------
// ALI SAFWAT

void DrawLineMidPoint(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int x = x1, y = y1;
    double dx = x2 - x1, dy = y2 - y1;
    SetPixel(hdc, x, y, c);
    if ((dx == 0 || dy / dx > 1) && dy > 0 && dx >= 0)
    {
        int d = 2 * dx - dy, d1 = 2 * dx, d2 = 2 * dx - 2 * dy;
        while (y != y2)
        {
            if (d <= 0)
            {
                y++;
                d += d1;
            }
            else
            {
                x++;
                y++;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy / dx >= 0 && dy / dx <= 1 && dy >= 0 && dx > 0)
    {
        int d = dx - 2 * dy, d1 = -2 * dy, d2 = 2 * dx - 2 * dy;
        while (x != x2)
        {
            if (d > 0)
            {
                x++;
                d += d1;
            }
            else
            {
                x++;
                y++;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy / dx < 0 && dy / dx >= -1 && dy <= 0 && dx>0)
    {
        int d = -dx - 2 * dy, d1 = -2 * dy, d2 = -2 * dx - 2 * dy;
        while (x != x2)
        {
            if (d <= 0)
            {
                x++;
                d += d1;
            }
            else
            {
                x++;
                y--;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx == 0 || dy / dx < -1) && dy < 0 && dx >= 0)
    {
        int d = -2 * dx - dy, d1 = -2 * dx, d2 = -2 * dx - 2 * dy;
        while (y != y2)
        {
            if (d > 0)
            {
                y--;
                d += d1;
            }
            else
            {
                x++;
                y--;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx == 0 || dy / dx > 1) && dy < 0 && dx <= 0)
    {
        int d = -2 * dx + dy, d1 = -2 * dx, d2 = -2 * dx + 2 * dy;
        while (y != y2)
        {
            if (d <= 0)
            {
                y--;
                d += d1;
            }
            else
            {
                x--;
                y--;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy / dx >= 0 && dy / dx <= 1 && dy <= 0 && dx < 0)
    {
        int d = -dx + 2 * dy, d1 = 2 * dy, d2 = -2 * dx + 2 * dy;
        while (x != x2)
        {
            if (d > 0)
            {
                x--;
                d += d1;
            }
            else
            {
                x--;
                y--;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if (dy / dx < 0 && dy / dx >= -1 && dy >= 0 && dx < 0)
    {
        int d = dx + 2 * dy, d1 = 2 * dy, d2 = 2 * dx + 2 * dy;
        while (x != x2)
        {
            if (d <= 0)
            {
                x--;
                d += d1;
            }
            else
            {
                x--;
                y++;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
    else if ((dx == 0 || dy / dx < -1) && dy > 0 && dx <= 0)
    {
        int d = 2 * dx + dy, d1 = 2 * dx, d2 = 2 * dx + 2 * dy;
        while (y != y2)
        {
            if (d > 0)
            {
                y++;
                d += d1;
            }
            else
            {
                x--;
                y++;
                d += d2;
            }
            SetPixel(hdc, x, y, c);
        }
    }
}
//struct POINT {
//    double x, y;
//};
void low_point(HDC hdc, int x0, int y0, int x1, int y1, COLORREF c) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int yi = 1;
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }
    int D = (2 * dy) - dx;
    int y = y0;

    for (int x = x0; x < x1; x++) {
        SetPixel(hdc, x, y, c);
        if (D > 0) {
            y = y + yi;
            D = D + (2 * (dy - dx));
        } else
            D = D + 2 * dy;
    }
}
void high_point(HDC hdc, int x0, int y0, int x1, int y1, COLORREF c) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int xi = 1;
    if (dx < 0) {
        xi = -1;
        dx = -dx;
    }
    int D = (2 * dx) - dy;
    int x = x0;

    for (int y = y0; y < y1; y++) {
        SetPixel(hdc, x, y, c);
        if (D > 0) {
            x = x + xi;
            D = D + (2 * (dx - dy));
        } else
            D = D + 2 * dx;
    }
}
void SetLine(HDC hdc, int x0, int y0, int x1, int y1, COLORREF c) {
    if (abs(y1 - y0) < abs(x1 - x0)) {
        if (x0 > x1)
            low_point(hdc, x1, y1, x0, y0, c);
        else
            high_point(hdc, x0, y0, x1, y1, c);
    } else {
        if (y0 > y1)
            low_point(hdc, x1, y1, x0, y0, c);
        else
            high_point(hdc, x0, y0, x1, y1, c);
    }
}
// point clipping
void PointClipping(HDC hdc, int x, int y, int xleft, int ytop, int xright, int ybottom, COLORREF color)
{
    if (x >= xleft && x <= xright && y >= ytop && y <= ybottom)
        SetPixel(hdc, x, y, color);
}
// line clipping
union OutCode {
    unsigned All: 4;
    struct {
        unsigned left: 1, top: 1, right: 1, bottom: 1;
    };
};
OutCode GetOutCode(POINT p, int xleft, int ytop, int xright, int ybottom) {
    OutCode out;
    out.All = 0;
    if (p.x < xleft)
        out.left = 1;
    else if (p.x > xright)
        out.right = 1;
    if (p.y > ytop)
        out.top = 1;
    else if (p.y < ybottom)
        out.bottom = 1;
    return out;
}
void vIntersection(POINT p1, POINT p2, int xEdge, POINT* P) {
    P->y = p1.y + (xEdge - p1.x) * (p2.y - p1.y) / (p2.x - p1.x);
    P->x = xEdge;
}
void HIntersection(POINT p1, POINT p2, int yEdge, POINT* P) {
    P->x = p1.x + (yEdge - p1.y) * (p2.x - p1.x) / (p2.y - p1.y);
    P->y = yEdge;
}

void LineClipping(HDC hdc, POINT p1, POINT p2, int xleft, int ytop, int xright, int ybottom) {
    OutCode out1 = GetOutCode(p1, xleft, ytop, xright, ybottom);
    OutCode out2 = GetOutCode(p2, xleft, ytop, xright, ybottom);

    while ((out1.All || out2.All) && !(out1.All & out2.All)) {
        if (out1.All) {
            if (out1.left)
                vIntersection(p1, p2, xleft, &p1);
            else if (out1.right)
                vIntersection(p1, p2, xright, &p1);
            else if (out1.top)
                HIntersection(p1, p2, ytop, &p1);
            else
                HIntersection(p1, p2, ybottom, &p1);
            out1 = GetOutCode(p1, xleft, ytop, xright, ybottom);
        }
        else {
            if (out2.left)
                vIntersection(p1, p2, xleft, &p2);
            else if (out2.right)
                vIntersection(p1, p2, xright, &p2);
            else if (out2.top)
                HIntersection(p1, p2, ytop, &p2);
            else
                HIntersection(p1, p2, ybottom, &p2);
            out2 = GetOutCode(p2, xleft, ytop, xright, ybottom);
        }
    }
    if (!out1.All && !out2.All) {
        SetLine(hdc, p1.x, p1.y, p2.x, p2.y, RGB(100, 100, 100));
    }
}

// polygon clipping
struct Vertex {
    double x, y;

    Vertex(int x1 = 0, int y1 = 0) {
        x = x1;
        y = y1;
    }
};
typedef vector<Vertex> VertexList;
typedef bool (*IsInFunc)(Vertex &v, int edge);
typedef Vertex (*IntersectFunc)(Vertex &v1, Vertex &v2, int edge);

VertexList ClipWithEdge(VertexList p, int edge, IsInFunc In, IntersectFunc Intersect) {
    VertexList OutList;
    Vertex v1 = p[p.size() - 1];
    bool v1_in = In(v1, edge);
    for (int i = 0; i < (int) p.size(); i++) {
        Vertex v2 = p[i];
        bool v2_in = In(v2, edge);
        if (!v1_in && v2_in) {

            OutList.push_back(Intersect(v1, v2, edge));
            OutList.push_back(v2);
        } else if (v1_in && v2_in) OutList.push_back(v2);
        else if (v1_in) OutList.push_back(Intersect(v1, v2, edge));
        v1 = v2;
        v1_in = v2_in;
    }
    return OutList;
}
bool InLeft(Vertex &v, int edge) {
    return v.x >= edge;
}
bool InRight(Vertex &v, int edge) {
    return v.x <= edge;
}
bool InTop(Vertex &v, int edge) {
    return v.y >= edge;
}
bool InBottom(Vertex &v, int edge) {
    return v.y <= edge;
}
Vertex VIntersect(Vertex &v1, Vertex &v2, int xedge) {
    Vertex res;
    res.x = xedge;
    res.y = v1.y + (xedge - v1.x) * (v2.y - v1.y) / (v2.x - v1.x);
    return res;
}
Vertex HIntersect(Vertex &v1, Vertex &v2, int yedge) {
    Vertex res;
    res.y = yedge;
    res.x = v1 .x + (yedge - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);
    return res;
}
void PolygonClip(HDC hdc, POINT *p, int n, int xleft, int ytop, int xright, int ybottom) {
    VertexList vlist;
    for (int i = 0; i < n; i++)
        vlist.push_back(Vertex(p[i].x, p[i].y));
    vlist = ClipWithEdge(vlist, xleft, InLeft, VIntersect);
    vlist = ClipWithEdge(vlist, ytop, InTop, HIntersect);
    vlist = ClipWithEdge(vlist, xright, InRight, VIntersect);
    vlist = ClipWithEdge(vlist, ybottom, InBottom, HIntersect);
    Vertex v1 = vlist[vlist.size() - 1];
    for (int i = 0; i < (int) vlist.size(); i++) {
        Vertex v2 = vlist[i];
        MoveToEx(hdc, Round(v1.x), Round(v1.y), NULL);
        LineTo(hdc, Round(v2.x), Round(v2.y));
        v1 = v2;
    }
}
void Square_clipping(HDC hdc, int top, int left, int R, COLORREF c)
{

    int right = left + R;
    int botton = top + R;

    if (botton < top) {
        swap(top, botton);
    }
    if (right < left) {
        swap(left, right);
    }

    DrawLineMidPoint(hdc, right, top, left, top, c);
    DrawLineMidPoint(hdc, right, botton, left, botton, c);
    DrawLineMidPoint(hdc, right, top, right, botton, c);
    DrawLineMidPoint(hdc, left, top, left, botton, c);
}
//------------------------------------------------------------------------------------------------------------------------------------------
// ABDO & AHMED
void swap(int &x, int &y) {
    int temp = x;
    x = y;
    y = temp;
}
int max(int x, int y) {
    return (x > y ? x : y);
}
int distance(int xc, int yc, int xr, int yr) {
    return (int) sqrt(pow(abs(xr - xc), 2) + pow(abs(yr - yc), 2));
}
void draw8Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c) {
    SetPixel(hdc, xc + x, yc + y, c);
    SetPixel(hdc, xc - x, yc + y, c);
    SetPixel(hdc, xc + x, yc - y, c);
    SetPixel(hdc, xc - x, yc - y, c);

    SetPixel(hdc, xc + y, yc + x, c);
    SetPixel(hdc, xc - y, yc + x, c);
    SetPixel(hdc, xc + y, yc - x, c);
    SetPixel(hdc, xc - y, yc - x, c);
}
void midPointCircle(HDC hdc, int xc, int yc, int rad, COLORREF color) {
    int x = 0;
    int y = rad;
    int d = 1 - rad;
    while (x <= y) {
        if (d <= 0) {
            x++;
            d += (2 * x) + 3;
        } else {
            x++;
            y--;
            d += 2 * (x - y) + 5;
        }
        draw8Points(hdc, xc, yc, x, y, color);
    }
}
void dda_line(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    // this algorithm use no * operations to reduce time consumed in the naive/direct algorithm
    // the idea is that : we can get the change and add it every time to the last co-ordinate (instead of calculating the co-ordinate)
    int dx = x2 - x1;
    int dy = y2 - y1;
    double m = 1.0 * dy / dx;   // slope
    double m2 = 1.0 / m;      // slope inverse
    if (abs(dx) >= abs(dy)) {
        // change in x > change in y
        if (x1 > x2) {
            // swap points
            swap(x1, x2);
            swap(y1, y2);
        }
        // x is independent , and y depends on x value
        // from the line equation : y = y1 + (x-x1)*m
        double yChange = m;    // newY - oldY  =  (y1+ (x+1 - x1)*m )  -  (y1 + (x-x1)*m) = m   ( change = (y at x+1) - (y at x) )

        // note that if y is int , the line will be drawn as a straight line ( because the value of m ( fraction ) will be ignored every time )
        double y = y1;
        for (int i = x1; i <= x2; ++i) {
            SetPixel(hdc, i, round(y), color);
            y += yChange;
        }
    } else {
        // change in x < change in y
        if (y1 > y2) {
            swap(y1, y2);
            swap(x1, x2);
        }
        // y is independent , x depends on y value
        // from line equation : x = x1 + (y-y1)*m2
        double xChange = m2;     // newX - oldX = (x1 + (y+1 - y1)*m2) - (x1 + (y - y1)*m2) = m2
        double x = x1;
        for (int i = y1; i <= y2; ++i) {
            SetPixel(hdc, round(x), i, color);
            x += xChange;
        }
    }
}
bool valid(int x, int y, int xc, int yc, int rad) {
    // if x inside the circle return true , otherwise return false ;
    int rightSide = rad * rad;

    int xLen = abs(x - xc);
    xLen *= xLen;

    int yLen = abs(y - yc);
    yLen *= yLen;
    int leftSide = (xLen + yLen);
    return (leftSide < rightSide);
}
void convexFillCircle(HDC hdc, int x1,int y1, int radius, int quarter, COLORREF color) {
    POINT center ;
    center.x = x1 ;
    center.y = y1 ;
    int boundary = center.y + radius;
    int yInc;
    int xInc;
    if (quarter == 1) {
        xInc = 1;
        yInc = -1;
    } else if (quarter == 2) {
        xInc = 1;
        yInc = 1;
    } else if (quarter == 3) {
        xInc = -1;
        yInc = 1;
    } else {
        xInc = -1;
        yInc = -1;
    }
    for (int y = center.y; y <= boundary; y += yInc) {
        for (int x = center.x; valid(x, y, center.x, center.y, radius); x += xInc)
            SetPixel(hdc, x, y, color);
    }
}
struct Vector2 {
    double x, y;

    Vector2(double a = 0, double b = 0) {
        x = a;
        y = b;
    }
};
class Vector4 {
    double v[4];
public:
    Vector4(double a = 0, double b = 0, double c = 0, double d = 0) {
        v[0] = a;
        v[1] = b;
        v[2] = c;
        v[3] = d;
    }

    Vector4(double a[]) {
        memcpy(v, a, 4 * sizeof(double));
    }

    double &operator[](int i) {
        return v[i];
    }
};
class Matrix4 {
    Vector4 M[4];
public:
    Matrix4(double A[]) {
        memcpy(M, A, 16 * sizeof(double));
    }

    Vector4 &operator[](int i) {
        return M[i];
    }
};
Vector4 operator*(Matrix4 M, Vector4 &b) {
    Vector4 res;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            res[i] += M[i][j] * b[j];

    return res;
}
double DotProduct(Vector4 &a, Vector4 &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}
Vector4 GetHermiteCoeff(double x0, double s0, double x1, double s1) {
    static double H[16] = {2, 1, -2, 1, -3, -2, 3, -1, 0, 1, 0, 0, 1, 0, 0, 0};
    static Matrix4 basis(H);
    Vector4 v(x0, s0, x1, s1);
    return basis * v;
}
void DrawHermiteCurve(HDC hdc, Vector2 &P0, Vector2 &T0, Vector2 &P1, Vector2 &T1, int numpoints, COLORREF color) {
    Vector4 xcoeff = GetHermiteCoeff(P0.x, T0.x, P1.x, T1.x);
    Vector4 ycoeff = GetHermiteCoeff(P0.y, T0.y, P1.y, T1.y);
    if (numpoints < 2)return;
    double dt = 1.0 / (numpoints - 1);
    for (double t = 0; t <= 1; t += dt) {
        Vector4 vt;
        vt[3] = 1;
        for (int i = 2; i >= 0; i--)vt[i] = vt[i + 1] * t;
        int x = round(DotProduct(xcoeff, vt));
        int y = round(DotProduct(ycoeff, vt));
        SetPixel(hdc, x, y, color);
    }
}
void DrawBezierCurve(HDC hdc, Vector2 &P0, Vector2 &P1, Vector2 &P2, Vector2 &P3, int numpoints, COLORREF color) {
    Vector2 T0(3 * (P1.x - P0.x), 3 * (P1.y - P0.y));
    Vector2 T1(3 * (P3.x - P2.x), 3 * (P3.y - P2.y));
    DrawHermiteCurve(hdc, P0, T0, P3, T1, numpoints, color);
}
void DrawCardinalSpline(HDC hdc, Vector2 P[], int n, double c, int numpix, COLORREF color) {
    double c1 = 1 - c;
    Vector2 T0(c1 * (P[2].x - P[0].x), c1 * (P[2].y - P[0].y));
    for (int i = 2; i < n - 1; i++) {
        Vector2 T1(c1 * (P[i + 1].x - P[i - 1].x), c1 * (P[i + 1].y - P[i - 1].y));
        DrawHermiteCurve(hdc, P[i - 1], T0, P[i], T1, numpix, color);
        T0 = T1;
    }
}
void fillWithBezier(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    if (x2 < x1) {
        int temp = x2;
        x2 = x1;
        x1 = temp;
    }
    if (y2 < y1) {
        int temp = y2;
        y2 = y1;
        y1 = temp;
    }
    for (int y = y1; y <= y2; y += 10) {
        for (int x = x1; x <= x2; x += 10) {
            Vector2 p0, t0, p1, t1;
            p0.x = x, p0.y = y;
            p1.x = x + 30, p1.y = y;
            t0.x = x + 30, t0.y = y + 30;
            t1 = t0;
            if (y + 30 <= y2 && x + 30 <= x2)
                DrawBezierCurve(hdc, p0, t0, p1, t1, 100, color);
        }
    }
}
void fillWithHermite(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    if (x2 < x1) {
        int temp = x2;
        x2 = x1;
        x1 = temp;
    }
    if (y2 < y1) {
        int temp = y2;
        y2 = y1;
        y1 = temp;
    }
    for (int y = y1 + 3; y <= y2; y += 3) {
        for (int x = x1; x <= x2; x += 3) {
            Vector2 p0, t0, p1, t1;
            p0.x = x, p0.y = y;
            t0.x = x + 2, t0.y = y + 2;
            p1.x = x + 5, p1.y = y + 3;
            t1.x = x + 2, t1.y = y;
            if (y + 12 <= y2 && x + 12 <= x2)
                DrawHermiteCurve(hdc, p0, t0, p1, t1, 100, color);
        }
    }
}
void square(HDC hdc, int x1, int y1, int len, COLORREF color, int border) {
    int x2=x1+len,y2=y1+len;
    int width = abs(x1 - x2);
    int height = abs(y1 - y2);
    int cnt = 0;
    if (x2 < x1) {
        int temp = x2;
        x2 = x1;
        x1 = temp;
    }
    if (y2 < y1) {
        int temp = y2;
        y2 = y1;
        y1 = temp;
    }
    for (int x = x1; x <= x1 + width; x++) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                x += 3;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                x += 5;
                continue;
            }
        }
        SetPixel(hdc, x, y1, color);
        SetPixel(hdc, x, y1 + height, color);
        cnt++;
    }
    cnt = 0;
    for (int y = y1; y <= y1 + height; y++) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                y += 3;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                y += 5;
                continue;
            }
        }
        SetPixel(hdc, x1, y, color);
        SetPixel(hdc, x1 + width, y, color);
        cnt++;
    }
}
void rrectangle(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color, int border) {
    int width = abs(x1 - x2);
    int height = abs(y1 - y2);
    int cnt = 0;
    if (x2 < x1) {
        int temp = x2;
        x2 = x1;
        x1 = temp;
    }
    if (y2 < y1) {
        int temp = y2;
        y2 = y1;
        y1 = temp;
    }
    for (int x = x1; x <= x1 + width; x++) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                x += 3;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                x += 5;
                continue;
            }
        }
        SetPixel(hdc, x, y1, color);
        SetPixel(hdc, x, y1 + height, color);
        cnt++;
    }
    cnt = 0;
    for (int y = y1; y <= y1 + height; y++) {
        if (border == 2) {  // dotted border
            if (cnt % 3 == 0) {
                cnt++;
                y += 3;
                continue;
            }
        } else if (border == 3) { // dashed border
            if (cnt % 10 == 0) {
                cnt++;
                y += 5;
                continue;
            }
        }
        SetPixel(hdc, x1, y, color);
        SetPixel(hdc, x1 + width, y, color);
        cnt++;
    }
}
HMENU hmenu;
void addmenu(HWND hwnd) {
    hmenu = CreateMenu();
    AppendMenu(hmenu, MF_STRING, 1, "save");
    AppendMenu(hmenu, MF_STRING, 2, "clear");
    AppendMenu(hmenu, MF_STRING, 3, "load");
    SetMenu(hwnd, hmenu);
}

bool saved = false;
//------------------------------------------------------------------------------------------------------------------------------------------
// ABDO
bool vis[1000][1000];
class point{
public:
    int x,y;
    point(int x=0,int y=0):x(x),y(y){}
};
void floodfill(HDC hdc, int x, int y, COLORREF bc, COLORREF fc) {
    if (vis[x][y])return;
    vis[x][y] = true;
    COLORREF c = GetPixel(hdc, x, y);
    if (c == bc || c == fc)return;
    SetPixel(hdc, x, y, fc);
    floodfill(hdc, x + 1, y, bc, fc);
    floodfill(hdc, x - 1, y, bc, fc);
    floodfill(hdc, x, y - 1, bc, fc);
    floodfill(hdc, x, y + 1, bc, fc);
}
void floodfill_nonrecursive(HDC hdc,int x,int y,COLORREF bc,COLORREF fc){
    stack<point>st;
    st.push(point(x,y));
    cout<<"hi"<<endl;
    while (!st.empty()){
        cout<<"inside"<<endl;
        point p=st.top();
        st.pop();
        COLORREF c= GetPixel(hdc,p.x,p.y);
        if(c==bc||c==fc)continue;
        SetPixel(hdc,x,y,fc);
    cout<<"puah"<<endl;
        st.push(point(p.x-1,p.y));
        st.push(point(p.x+1,p.y));
        st.push(point(p.x,p.y-1));
        st.push(point(p.x,p.y+1));
    }
}
void drawlineDDA(HDC hdc, int x1, int y1, int x2, int y2, COLORREF color) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    SetPixel(hdc, x1, y1, color);
    if (abs(dx) >= abs(dy)) {
        int x = x1, change1;
        if (dx > 0) {
            change1 = 1;
        } else {
            change1 = -1;
        }
        double y = y1, change2 = (double) dy / dx * change1;
        while (x != x2) {
            x += change1;
            y += change2;
            SetPixel(hdc, x, Round(y), color);
        }
    } else {
        int y = y1, change2;
        if (dy > 0) {
            change2 = 1;
        } else {
            change2 = -1;
        }
        double x = x1, change1 = (double) dx / dy * change2;
        while (y != y2) {
            x += change1;
            y += change2;
            SetPixel(hdc, Round(x), y, color);
        }
    }
}
void midPointCircle(HDC hdc, int xc, int yc, int radius, int q, COLORREF color) {
    int x = 0, y = radius;
    int d = 1 - radius;
    draw8Points(hdc, xc, yc, x, y, color);
    while (x <= y) {
        if (d <= 0) {
            d += (2 * x) + 3;
            x++;
        } else {
            d += 2 * (x - y) + 5;
            x++;
            y--;
        }
        draw8Points(hdc, xc, yc, x, y, color);
    }
    int x1, y1, x2, y2;
    //V
    x1 = xc, y1 = yc - radius, x2 = xc, y2 = yc + radius;
    drawlineDDA(hdc, x1, y1, x2, y2, color);
    //H
    x1 = xc - radius, y1 = yc, x2 = xc + radius, y2 = yc;
    drawlineDDA(hdc, x1, y1, x2, y2, color);
}
void fill_quarter(HDC hdc, int xc, int yc,int radius, int q, COLORREF color) {
   // int x1,y1,x2,y2;
    //V
    //x1=xc, y1=yc - radius, x2=xc, y2=yc + radius;
    //drawlineDDA(hdc, x1, y1, x2, y2, color);
    //H
    //x1=xc - radius, y1=yc, x2=xc + radius, y2=yc;
    //drawlineDDA(hdc, x1, y1, x2 , y2, color);
    COLORREF color_inside = RGB(133, 87, 35);
    if (q == 1)
        floodfill(hdc, xc - 1, yc - 1, color, color_inside);
    else if (q == 2)
        floodfill(hdc, xc + 1, yc - 1, color, color_inside);
    else if (q == 3)
        floodfill(hdc, xc + 1, yc + 1, color, color_inside);
    else
        floodfill(hdc, xc - 1, yc + 1, color, color_inside);
}
void LineParametric(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
    int dx = x2 - x1;
    int dy = y2 - y1;
    int steps = max(abs(dx), abs(dy));

    float increment = 1.0 / steps;

    float x = x1;
    float y = y1;

    for (int i = 0; i <= steps; i++) {
        SetPixel(hdc, x, y, c);
        x += dx * increment;
        y += dy * increment;
    }
}
vector<POINT> points;
void CircleFasterBresenham(HDC hdc, int xc, int yc, int R, COLORREF color) {
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y) {
		if (d < 0) {
			d += c1;
			c2 += 2;
		}
		else {

			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
		points.push_back({ x, y });
	}
}
void FilledCircleQuarterWithCircles(HDC hdc, int xc, int yc, int radius, int quarter, COLORREF c) {
    // Calculate the squared radius
    int radiusSquared = radius * radius;
    CircleFasterBresenham(hdc, xc, yc, radius, c);

    // Vector to store the center points of the smaller circles
    std::vector<POINT> centerPoints;
    int startX, startY, endX, endY;
    switch (quarter) {
    case 1:  // Upper-right quarter
        startX = 0;
        startY = -radius;
        endX = radius;
        endY = 0;
        break;
    case 2:  // Upper-left quarter
        startX = -radius;
        startY = -radius;
        endX = 0;
        endY = 0;
        break;
    case 3:  // Lower-left quarter
        startX = -radius;
        startY = 0;
        endX = 0;
        endY = radius;
        break;
    case 4:  // Lower-right quarter
        startX = 0;
        startY = 0;
        endX = radius;
        endY = radius;
        break;
    default:
        std::cout << "Invalid quarter specified!" << std::endl;
        return;
    }

    // Iterate over the bounding box of the filling quarter
    for (int y = startY; y <= endY; y++) {
        for (int x = startX; x <= endX; x++) {
            if (x * x + y * y <= radiusSquared) {
                centerPoints.push_back({ xc + x, yc + y });
            }
        }
    }

    // Draw smaller circles inside the filled region
    int smallerRadius = radius / 10;
    for (const POINT& center : centerPoints) {
        CircleFasterBresenham(hdc, center.x, center.y, smallerRadius, c);
    }
}



/*--------------Start Filling-----------------------*/

//------------------------- Line clipping --------------------------//


/*--------------------save and load----------------------*/


/*---------------------------------------------------*/
/*--------------------------spiling-------------------*/

/*-----------------------------------------------*/
/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure(HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow) {
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof(WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon(NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor(NULL, IDC_CROSS);///////////////MOUCE
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = (HBRUSH) RGB(255, 255, 0);

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx(&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx(
            0,                   /* Extended possibilites for variation */
            szClassName,         /* Classname */
            _T("Project Graphics"),       /* Title Text */
            WS_OVERLAPPEDWINDOW, /* default window */
            CW_USEDEFAULT,       /* Windows decides the position */
            CW_USEDEFAULT,       /* where the window ends up on the screen */
            544,                 /* The programs width */
            375,                 /* and height in pixels */
            HWND_DESKTOP,        /* The window is a child-window to desktop */
            NULL,                /* No menu */
            hThisInstance,       /* Program Instance handler */
            NULL                 /* No Window Creation data */
    );

    /* Make the window visible on the screen */
    ShowWindow(hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage(&messages, NULL, 0, 0)) {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}


/*  This function is called by the Windows function DispatchMessage()  */
COLORREF color = RGB(0, 0, 0);
int counter = 0;
static int floodCnt=0 ;
Vector2 P[7];
POINT p1, p2;
int case_number = 0, quarter = 0;
int xleft, ytop, xright, ybottom;
int index = 0;
static int fillCounter = 0 ;
static int border;
static int circleCnt=0 ;
LRESULT CALLBACK WindowProcedure(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    HDC hdc = GetDC(hwnd);
    static int x1, x2, x3 = 0;
    static int y1, y2, y3, x4 = 0, y4 = 0;
    int r, r2 = 0;
    static int R = 0;
    switch (message)                  /* handle the messages */
    {
        case WM_LBUTTONDOWN:
                if(counter==0){
                    x1 = LOWORD(lParam);
                    y1 = HIWORD(lParam);
                    counter++ ;
                }
                else if(counter==1){
                    x2 = LOWORD(lParam);
                    y2 = HIWORD(lParam);
                    counter++ ;
                }
                else{
                    CircleFasterBresenham(hdc,x1,y1,distance(x1,y1,x2,y2),color);
                }


        case WM_DESTROY:
            PostQuitMessage(0);       /* send a WM_QUIT to the message queue */
            break;

        default:                      /* for messages that we don't deal with */
            return DefWindowProc(hwnd, message, wParam, lParam);
    }

    return 0;
}
