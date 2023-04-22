// ACM_ICPC算法与实现(计算几何专题)立体计算几何
#include <bits/stdc++.h>
using namespace std;

const double eps = 1e-8;
const int N = 100010;
const double PI = acos(-1);

inline int cmp(double a)
{
    return a < -eps ? -1 : a > eps;
}

inline double Sqr(double a)
{
    return a * a;
}

inline double Sqrt(double a)
{
    return a <= 0 ? 0 : sqrt(a);
}

class Point_3
{
public:
    double x, y, z;
    Point_3() {}
    Point_3(double x, double y, double z) : x(x), y(y), z(z) {}

    // 向量长度
    double Length() const
    {
        return Sqrt(Sqr(x) + Sqr(y) + Sqr(z));
    }

    Point_3 Unit() const;

    bool operator<(const Point_3 &Range) const
    {
        if (x == Range.x)
        {
            if (y == Range.y)
                return z < Range.z;
            else
                return y < Range.y;
        }
        else
            return x < Range.x;
    }
};

Point_3 operator+(const Point_3 &a, const Point_3 &b)
{
    return Point_3(a.x + b.x, a.y + b.y, a.z + b.z);
}
Point_3 operator-(const Point_3 &a, const Point_3 &b)
{
    return Point_3(a.x - b.x, a.y - b.y, a.z - b.z);
}
Point_3 operator*(const Point_3 &a, const Point_3 &b)
{
    return Point_3(a.y * b.z - b.y * a.z, -a.x * b.z + a.z * b.x, a.x * b.y - a.y * b.x);
}
Point_3 operator*(double t, const Point_3 &a)
{
    return Point_3(t * a.x, t * a.y, t * a.z);
}
Point_3 operator*(const Point_3 &a, double t)
{
    return Point_3(t * a.x, t * a.y, t * a.z);
}
Point_3 operator/(const Point_3 &a, double t)
{
    return Point_3(a.x / t, a.y / t, a.z / t);
}

// 返回单位化的向量
Point_3 Point_3::Unit() const
{
    return *this / Length();
}

// 向量a和b的点积
double Dot(const Point_3 &a, const Point_3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// 向量a,b,c的混合积.返回值除以6就是四面体体积
double Mix(const Point_3 &a, const Point_3 &b, const Point_3 &c)
{
    return Dot(a, b * c);
}

// 三维两点距离
double distance(const Point_3 &a, const Point_3 &b)
{
    return Sqrt(Sqr(a.x - b.x) + Sqr(a.y - b.y) + Sqr(a.z - b.z));
}

// 三维直线类
struct Line_3
{
    Point_3 a, b;
    Line_3() {}
    Line_3(Point_3 a, Point_3 b) : a(a), b(b) {}
};

// 线段长度
double vlen(Point_3 P)
{
    return P.Length();
}
// 零值函数
bool zero(double x)
{
    return fabs(x) < eps;
}
// 判断三点共线
int dots_inline(Point_3 p1, Point_3 p2, Point_3 p3)
{
    return vlen((p2 - p1) * (p2 - p3)) < eps;
}
// 判断点在线段内(包含端点)
int dot_online_in(Point_3 p, Line_3 l)
{
    return zero(vlen((p - l.a) * (p - l.b))) && (l.a.x - p.x) * (l.b.x - p.x) < eps && (l.a.y - p.y) * (l.b.y - p.y) < eps && (l.a.z - p.z) * (l.b.z - p.z) < eps;
}
// 判断点在线段内(不包含端点)
int dot_online_ex(Point_3 p, Line_3 l)
{
    return dot_online_in(p, l) && (!zero(p.x - l.a.x) || !zero(p.y - l.a.y) || !zero(p.z - l.a.z)) && (!zero(p.x - l.b.x) || !zero(p.y - l.b.y) || !zero(p.z - l.b.z));
}
// 判断平面内两点在直线同侧
int same_side(Point_3 p1, Point_3 p2, Line_3 l)
{
    return Dot((l.a - l.b) * (p1 - l.b), (l.a - l.b) * (p2 - l.b)) > eps;
}
// 判断平面内两点在直线异侧
int opposite_side(Point_3 p1, Point_3 p2, Line_3 l)
{
    return Dot((l.a - l.b) * (p1 - l.b), (l.a - l.b) * (p2 - l.b)) < -eps;
}
// 判断两直线平行
int parallel(Line_3 u, Line_3 v)
{
    return vlen((u.a - u.b) * (v.a - v.b)) < eps;
}
// 判断两直线垂直
int perpendicular(Line_3 u, Line_3 v)
{
    return zero(Dot(u.a - u.b, v.a - v.b));
}

// 三维平面类
struct Plane_3
{
    Point_3 a, b, c;
    Plane_3() {}
    Plane_3(Point_3 a, Point_3 b, Point_3 c) : a(a), b(b), c(c) {}
};
// 平面法向量
Point_3 pvec(Point_3 s1, Point_3 s2, Point_3 s3)
{
    return (s1 - s2) * (s2 - s3);
}
// 判断四点共面
int dots_onplane(Point_3 a, Point_3 b, Point_3 c, Point_3 d)
{
    return zero(Dot(pvec(a, b, c), d - a));
}
// 直线版重载
int dots_onplane(Line_3 u, Line_3 v)
{
    return dots_onplane(u.a, u.b, v.a, v.b);
}
// 判断两条线断是否有交点(包含端点)
int intersect_in(Line_3 u, Line_3 v)
{
    if (!dots_onplane(u, v))
        return 0;
    if (!dots_inline(u.a, u.b, v.a) || !dots_inline(u.a, u.b, v.b))
        return !same_side(u.a, u.b, v) && !same_side(v.a, v.b, u);
    return dot_online_in(u.a, v) || dot_online_in(u.b, v) || dot_online_in(v.a, u) || dot_online_in(v.b, u);
}
// 判断两条线断是否有交点(不包含端点)
int intersect_ex(Line_3 u, Line_3 v)
{
    return dots_onplane(u, v) && opposite_side(u.a, u.b, v) && opposite_side(v.a, v.b, u);
}
// 求两直线交点(保证共面且不平行)
Point_3 intersection(Line_3 u, Line_3 v)
{
    Point_3 ret = u.a;
    double t = ((u.a.x - v.a.x) * (v.a.y - v.b.y) - (u.a.y - v.a.y) * (v.a.x - v.b.x)) / ((u.a.x - u.b.x) * (v.a.y - v.b.y) - (u.a.y - u.b.y) * (v.a.x - v.b.x));
    ret = ret + (u.b - u.a) * t;
    return ret;
}
// 点到直线距离
double ptoline(Point_3 p, Line_3 l)
{
    return vlen((p - l.a) * (l.b - l.a)) / distance(l.a, l.b);
}
// 直线到直线距离,平行时需特殊处理
double linetoline(Line_3 u, Line_3 v)
{
    Point_3 n = (u.a - u.b) * (v.a - v.b);
    return fabs(Dot(u.a - v.a, n) / vlen(n));
}
// 求两直线夹角的cos值
double angle_cos(Line_3 u, Line_3 v)
{
    return Dot(u.a - u.b, v.a - v.b) / vlen(u.a - u.b) / vlen(v.a - v.b);
}
// 判断一个点是否在三角形内(包含边界)
int dot_inplane_in(Point_3 p, Plane_3 s)
{
    return zero(vlen((s.a - s.b) * (s.a - s.c)) - vlen(((p - s.a) * (p - s.b))) - vlen((p - s.b) * (p - s.c)) - vlen((p - s.c) * (p - s.a)));
}
// 判断一个点是否在三角形内(不包含边界)
int dot_inplane_ex(Point_3 p, Plane_3 s)
{
    return dot_inplane_in(p, s) && vlen((p - s.a) * (p - s.b)) > eps && vlen((p - s.b) * (p - s.c)) > eps && vlen((p - s.c) * (p - s.a)) > eps;
}
void solve()
{
}
int n;
Point_3 point[N];

int main()
{
    // ios::sync_with_stdio(false);
    // cin.tie(0) , cout.tie(0);
    cin >> n;
    for (int i = 0; i < n; i++)
        cin >> point[i].x >> point[i].y >> point[i].z;

    solve();

    system("pause");
    return 0;
}