// ACM_ICPC算法与实现(计算几何专题)平面计算几何
#include <bits/stdc++.h>
#include <bits/extc++.h>
using namespace std;
using namespace __gnu_pbds; // 红黑树、平衡树、字符串哈希标准库
// 3.1多边形---
const double eps = 1e-8;

// 符号函数
int cmp(double x)
{
    if (fabs(x) < eps)
        return 0;
    if (x < 0)
        return -1;
    return 1;
}

const double pi = acos(-1.0);

// 计算一个数的平方
inline double sqr(double x)
{
    return x * x;
}

// 二维点类
struct point
{
    double x, y;
    point() {}
    point(double a, double b) : x(a), y(b) {}
    // 读入数据
    void input()
    {
        scanf("%lf%lf", &x, &y);
    }
    // 运算符重载
    friend point operator+(const point &a, const point &b)
    {
        return point(a.x + b.x, a.y + b.y);
    }
    friend point operator-(const point &a, const point &b)
    {
        return point(a.x - b.x, a.y - b.y);
    }
    friend bool operator==(const point &a, const point &b)
    {
        return cmp(a.x - b.x) == 0 && cmp(a.y - b.y) == 0;
    }
    friend point operator*(const point &a, const double &b)
    {
        return point(a.x * b, a.y * b);
    }
    friend point operator*(const double &a, const point &b)
    {
        return point(a * b.x, a * b.y);
    }
    friend point operator/(const point &a, const double &b)
    {
        return point(a.x / b, a.y / b);
    }
    double norm()
    {
        return sqrt(sqr(x) + sqr(y));
    }
};

// 计算两个向量的叉积
double det(const point &a, const point &b)
{
    return a.x * b.y - a.y * b.x;
}
// 计算两个向量的点积
double dot(const point &a, const point &b)
{
    return a.x * b.x + a.y * b.y;
}
// 计算两个点的距离
double dist(const point &a, const point &b)
{
    return (a - b).norm();
}
// op向量绕原点逆时针旋转A(弧度)
point rotate_point(const point &p, double A)
{
    double tx = p.x, ty = p.y;
    return point(tx * cos(A) - ty * sin(A), tx * sin(A) + ty * cos(A));
}

// 线段类
struct line
{
    point a, b;
    line() {}
    line(point x, point y) : a(x), b(y) {}
};

// 用a,b两个点生成一个线段或直线
line point_make_line(const point a, const point b)
{
    return line(a, b);
}
// 求P点到线段st的距离
double dis_point_segment(const point p, const point s, const point t)
{
    if (cmp(dot(p - s, t - s)) < 0)
        return (p - s).norm();
    if (cmp(dot(p - t, s - t)) < 0)
        return (p - t).norm();
    return fabs(det(s - p, t - p) / dist(s, t));
}
// 求P点到线段st的垂足,保存在cp中
void PointProjLine(const point p, const point s, const point t, point &cp)
{
    double r = dot((t - s), (p - s)) / dot(t - s, t - s);
    cp = s + r * (t - s);
}
// 判断P点是否在线段st上(包括端点)
bool PointOnSegment(point p, point s, point t)
{
    return cmp(det(p - s, t - s)) == 0 && cmp(dot(p - s, p - t)) <= 0;
}
// 判断a和b是否平行
bool parallel(line a, line b)
{
    return !cmp(det(a.a - a.b, b.a - b.b));
}
// 判断a和b是否相交，如果相交则返回true并且交点保存在res中
bool line_make_point(line a, line b, point &res)
{
    if (parallel(a, b))
        return false;
    double s1 = det(a.a - b.a, b.b - b.a);
    double s2 = det(a.b - b.a, b.b - b.a);
    res = (s1 * a.b - s2 * a.a) / (s1 - s2);
    return true;
}
// 将直线a沿法向量方向平移距离len得到的直线f'r
line move_d(line a, const double &len)
{
    point d = a.b - a.a;
    d = d / d.norm();
    d = rotate_point(d, pi / 2);
    return line(a.a + d * len, a.b + d * len);
}

// 多边形类
const int maxn = 1010;
struct polygon
{
    // 成员变量
    int n;         // 多边形点数
    point a[maxn]; // 顺时针多边形顶点坐标
    polygon() {}

    // 读入数据
    void input()
    {
        scanf("%d", &n);
        for (int i = 0; i < n; i++)
            a[i].input();
    }

    // 成员函数

    // 计算多边形周长
    double perimeter()
    {
        double sum = 0;
        a[n] = a[0];
        for (int i = 0; i < n; i++)
            sum += (a[i + 1] - a[i]).norm();
        return sum;
    }
    // 计算多边形面积
    double area()
    {
        double sum = 0;
        a[n] = a[0];
        for (int i = 0; i < n; i++)
            sum += det(a[i + 1], a[i]);
        return sum / 2.;
    }
    // 判断点是否在多边形内部(0外1内2边界上)
    int Point_In(point t)
    {
        int num = 0;
        a[n] = a[0];
        for (int i = 0; i < n; i++)
        {
            if (PointOnSegment(t, a[i], a[i + 1]))
                return 2;
            int k = cmp(det(a[i + 1] - a[i], t - a[i]));
            int d1 = cmp(a[i].y - t.y);
            int d2 = cmp(a[i + 1].y - t.y);
            if (k > 0 && d1 <= 0 && d2 > 0)
                num++;
            if (k < 0 && d2 <= 0 && d1 > 0)
                num--;
        }
        return num != 0;
    }
    // 多边形的重心
    point MassCenter()
    {
        point ans = point(0, 0);
        if (cmp(area()) == 0)
            return ans;
        a[n] = a[0];
        for (int i = 0; i < n; i++)
            ans = ans + (a[i] + a[i + 1]) * det(a[i + 1], a[i]);
        return ans / area() / 6.;
    }
    // 多边形内格点数(Pick公式)
    int Border_Int_Point_Num()
    {
        int num = 0;
        a[n] = a[0];
        for (int i = 0; i < n; i++)
            num += __gcd(abs(int(a[i + 1].x - a[i].x)), abs(int(a[i + 1].y - a[i].y)));
        return num;
    }
    int Inside_Int_Point_Num()
    {
        return int(area()) + 1 - Border_Int_Point_Num() / 2;
    }
};

// 凸多边形类
struct polygon_convex
{
    vector<point> P;
    polygon_convex(int Size = 0)
    {
        P.resize(Size);
    }
};

bool comp_less(const point &a, const point &b)
{
    return cmp(a.x - b.x) < 0 || cmp(a.x - b.x) == 0 && cmp(a.y - b.y) < 0;
}
// 求凸包
polygon_convex convex_hull(vector<point> a)
{
    polygon_convex res(2 * a.size() + 5);
    sort(a.begin(), a.end(), comp_less);
    a.erase(unique(a.begin(), a.end()), a.end());
    int m = 0;
    for (int i = 0; i < a.size(); i++)
    {
        while (m > 1 && cmp(det(res.P[m - 1] - res.P[m - 2], a[i] - res.P[m - 2])) <= 0)
            m--;
        res.P[m++] = a[i];
    }
    int k = m;
    for (int i = int(a.size()) - 2; i >= 0; i--)
    {
        while (m > k && cmp(det(res.P[m - 1] - res.P[m - 2], a[i] - res.P[m - 2])) <= 0)
            m--;
        res.P[m++] = a[i];
    }
    res.P.resize(m);
    if (a.size() > 1)
        res.P.resize(m - 1);
    return res;
}
bool containOn(const polygon_convex &a, const point &b)
{
    int n = a.P.size();
#define next(i) ((i + 1) % n)
    int sign = 0;
    for (int i = 0; i < n; i++)
    {
        int x = cmp(det(a.P[i] - b, a.P[next(i)] - b));
        if (x)
        {
            if (sign)
            {
                if (sign != x)
                    return false;
            }
            else
                sign = x;
        }
    }
    return true;
}
// b是否在a的内部或边界上
int containOlogn(const polygon_convex &a, const point &b)
{
    int n = a.P.size();
    // 找到凸包内部的点g
    point g = (a.P[0] + a.P[n / 3] + a.P[2 * n / 3]) / 3.0;
    int l = 0, r = n;
    // 二分凸包 g - a.P[a] - a.P[b]
    while (l + 1 < r)
    {
        int mid = l + r >> 1;
        if (cmp(det(a.P[l] - g, a.P[mid] - g)) > 0)
        {
            if (cmp(det(a.P[l] - g, b - g)) >= 0 && cmp(det(a.P[mid] - g, b - g)) < 0)
                r = mid;
            else
                l = mid;
        }
        else
        {
            if (cmp(det(a.P[l] - g, b - g)) < 0 && cmp(det(a.P[mid] - g, b - g)) >= 0)
                l = mid;
            else
                r = mid;
        }
    }
    r %= n;
    int z = cmp(det(a.P[l] - g, b - g)) - 1;
    if (z == -2)
        return 1;
    return z;
}
// 凸多边形的直径(旋转卡壳)
double convex_diameter(polygon_convex &a, int &First, int &Second)
{
    vector<point> &p = a.P;
    int n = p.size();
    double maxd = 0.0;
    if (n == 1)
    {
        First = Second = 0;
        return maxd;
    }
#define next(i) ((i + 1) % n)
    for (int i = 0, j = 1; i < n; i++)
    {
        while (cmp(det(p[next(i)] - p[i], p[j] - p[i]) - det(p[next(i)] - p[i], p[next(j)] - p[i])) < 0)
            j = next(j);
        double d = dist(p[i], p[j]);
        if (d > maxd)
        {
            maxd = d;
            First = i, Second = j;
        }
        d = dist(p[next(i)], p[next(j)]);
        if (d > maxd)
        {
            maxd = d;
            First = i, Second = j;
        }
    }
    return maxd;
}

// 半平面类
struct halfPlane
{
    // ax + by + c <= 0
    double a, b, c;

    halfPlane(point p, point q)
    {
        a = p.y - q.y;
        b = q.x - p.x;
        c = det(p, q);
    }
    halfPlane(double aa, double bb, double cc)
    {
        a = aa, b = bb, c = cc;
    }
};

// 计算点a带入到直线方程中的函数值
double calc(halfPlane &L, point &a)
{
    return a.x * L.a + a.y * L.b + L.c;
}

// 求点a和b连线与半平面L的交点
point Intersect(point &a, point &b, halfPlane &L)
{
    point res;
    double t1 = calc(L, a), t2 = calc(L, b);
    res.x = (t2 * a.x - t1 * b.x) / (t2 - t1);
    res.y = (t2 * a.y - t1 * b.y) / (t2 - t1);
    return res;
}

// 半平面切割多边形
polygon_convex cut(polygon_convex &a, halfPlane &L)
{
    int n = a.P.size();
    polygon_convex res;
    for (int i = 0; i < n; i++)
    {
        if (calc(L, a.P[i]) < -eps)
            res.P.push_back(a.P[i]);
        else
        {
            int j = i - 1;
            if (j < 0)
                j = n - 1;
            if (calc(L, a.P[j]) < -eps)
                res.P.push_back(Intersect(a.P[j], a.P[i], L));
            j = i + 1;
            if (j == n)
                j = 0;
            if (calc(L, a.P[j]) < -eps)
                res.P.push_back(Intersect(a.P[i], a.P[j], L));
        }
    }
    return res;
}

// 半平面交
typedef complex<double> Point;
typedef pair<Point, Point> Halfplane;
const double EPS = 1e-12;
const double INF = 10000;

inline int sgn(double n)
{
    return fabs(n) < EPS ? 0 : (n < 0 ? -1 : 1);
}
inline double cross(Point a, Point b)
{
    return (conj(a) * b).imag();
}
inline double dot(Point a, Point b)
{
    return (conj(a) * b).real();
}
inline double satisfy(Point a, Halfplane p)
{
    return sgn(cross(a - p.first, p.second - p.first)) <= 0;
}
Point crosspoint(const Halfplane &a, const Halfplane &b)
{
    double k = cross(b.first - b.second, a.first - b.second);
    k = k / (k - cross(b.first - b.second, a.second - b.second));
    return a.first + (a.second - a.first) * k;
}
bool comp_Halfplane(const Halfplane &a, const Halfplane &b)
{
    int res = sgn(arg(a.second - a.first) - arg(b.second - b.first));
    return res == 0 ? satisfy(a.first, b) : res < 0;
}

vector<Point> halfplaneIntersection(vector<Halfplane> v)
{
    sort(v.begin(), v.end(), comp_Halfplane);
    deque<Halfplane> q;
    deque<Point> ans;
    q.push_back(v[0]);
    for (int i = 1; i < int(v.size()); i++)
    {
        if (sgn(arg(v[i].second - v[i].first) - arg(v[i - 1].second - v[i - 1].first)) == 0)
            continue;
        while (ans.size() > 0 && !satisfy(ans.back(), v[i]))
        {
            ans.pop_back();
            q.pop_back();
        }
        while (ans.size() > 0 && !satisfy(ans.front(), v[i]))
        {
            ans.pop_front();
            q.pop_front();
        }
        ans.push_back(crosspoint(q.back(), v[i]));
        q.push_back(v[i]);
    }
    while (ans.size() > 0 && !satisfy(ans.back(), q.front()))
    {
        ans.pop_back();
        q.pop_back();
    }
    while (ans.size() > 0 && !satisfy(ans.front(), q.back()))
    {
        ans.pop_front();
        q.pop_front();
    }
    ans.push_back(crosspoint(q.back(), q.front()));
    return vector<Point>(ans.begin(), ans.end());
}

// 凸多边形交
typedef vector<Point> Convex;
Convex convexIntersection(Convex v1, Convex v2)
{
    vector<Halfplane> h;
    for (int i = 0; i < int(v1.size()); i++)
        h.push_back(Halfplane(v1[i], v1[(i + 1) % v1.size()]));
    for (int i = 0; i < int(v2.size()); i++)
        h.push_back(Halfplane(v2[i], v2[(i + 1) % v2.size()]));
    return halfplaneIntersection(h);
}

// 多边形的核
const double inf = 10000;
polygon_convex core(polygon &a)
{
    polygon_convex res;
    res.P.push_back(point(-inf, -inf));
    res.P.push_back(point(inf, -inf));
    res.P.push_back(point(inf, inf));
    res.P.push_back(point(-inf, inf));

    int n = a.n;
    for (int i = 0; i < n; i++)
    {
        halfPlane L(a.a[i], a.a[(i + 1) % n]);
        res = cut(res, L);
    }
    return res;
}

// 凸多边形与直线交集
inline bool operator<(const point &a, const point &b)
{
    return a.y + eps < b.y || fabs(a.y - b.y) < eps && a.x + eps < b.x;
}
inline double getA(const point &a) // 凸包顺序下对应的角度也要递增
{
    double res = atan2(a.y, a.x);
    if (res < 0)
        res += 2 * pi;
    return res;
}

point p[maxn], hull[maxn];
int n;
double w[maxn], sum[maxn];

inline void GetHull() // 预处理一个逆时针的凸包
{
    sort(p + 1, p + n + 1);
    int N = 0;
    hull[++N] = p[1];
    for (int i = 2; i <= n; i++)
    {
        while (N > 1 && det(hull[N] - hull[N - 1], p[i] - hull[N - 1]) <= 0)
            N--;
        hull[++N] = p[i];
    }
    int bak = N;
    for (int i = n - 1; i >= 1; i--)
    {
        while (N > bak && det(hull[N] - hull[N - 1], p[i] - hull[N - 1]) <= 0)
            N--;
        hull[++N] = p[i];
    }
    n = N - 1;
    for (int i = 1; i <= n; i++)
        p[i + n] = p[i] = hull[i];
    p[n + n + 1] = p[1];
    for (int i = 1; i <= n; i++)
        w[i + n] = w[i] = getA(p[i + 1] - p[i]);
    sum[0] = 0;
    for (int i = 1; i <= 2 * n; i++)
        sum[i] = sum[i - 1] + det(p[i], p[i + 1]); // 预处理有向面积的前缀和
}

inline int Find(double x) // 找第一个角度>x的边
{
    if (x <= w[1] || x >= w[n])
        return 1;
    return (upper_bound(w + 1, w + n + 1, x) - (w + 1)) + 1;
}

point P, Q;
inline int getInter(int l, int r) // 找到第一个和p[1]不在P->Q向量同侧的点
{
    int sign;
    if (det(Q - P, p[1] - P) < 0)
        sign = -1;
    else
        sign = 1;
    while (l + 1 < r)
    {
        int mid = l + r >> 1;
        if (det(Q - P, p[mid] - P) * sign > 0)
            l = mid;
        else
            r = mid;
    }
    return r;
}

inline point Intersect(const point &a, const point &b, const point &c, const point &d) // 两直线求交点
{
    double s1 = det(c - a, b - a);
    double s2 = det(d - a, b - a);
    return (c * s2 - d * s1) / (s2 - s1);
}
inline bool solve(point P, point Q)
{
    int i = Find(getA(Q - P));
    int j = Find(getA(P - Q));
    // 两侧各找一点
    if (det(Q - P, p[i] - P) * det(Q - P, p[j] - P) >= 0) // 无交点
        return false;
    else
        return true;
}

// 3.2圆

// 精度符号sign函数
int dcmp(double k)
{
    return k < -EPS ? -1 : k > EPS ? 1
                                   : 0;
}
// 求根号(防止开负数根的情况)
double mysqrt(double n)
{
    return sqrt(max(0.0, n));
}
// 圆与直线求交
void circle_cross_line(point a, point b, point o, double r, point ret[], int &num)
{
    double x0 = o.x, y0 = o.y;
    double x1 = a.x, y1 = a.y;
    double x2 = b.x, y2 = b.y;
    double dx = x2 - x1, dy = y2 - y1;
    double A = dx * dx + dy * dy;
    double B = 2 * dx * (x1 - x0) + 2 * dy * (y1 - y0);
    double C = sqr(x1 - x0) + sqr(y1 - y0) - sqr(r);
    double delta = B * B - 4 * A * C;
    num = 0;
    if (dcmp(delta) >= 0)
    {
        double t1 = (-B - mysqrt(delta)) / (2 * A);
        double t2 = (-B + mysqrt(delta)) / (2 * A);
        ret[num++] = point(x1 + t1 * dx, y1 + t1 * dy);
        ret[num++] = point(x1 + t2 * dx, y1 + t2 * dy);
    }
}

// 圆与线段求交
void circle_cross_segment(point a, point b, point o, double r, point ret[], int &num)
{
    double x0 = o.x, y0 = o.y;
    double x1 = a.x, y1 = a.y;
    double x2 = b.x, y2 = b.y;
    double dx = x2 - x1, dy = y2 - y1;
    double A = dx * dx + dy * dy;
    double B = 2 * dx * (x1 - x0) + 2 * dy * (y1 - y0);
    double C = sqr(x1 - x0) + sqr(y1 - y0) - sqr(r);
    double delta = B * B - 4 * A * C;
    num = 0;
    if (dcmp(delta) >= 0)
    {
        double t1 = (-B - mysqrt(delta)) / (2 * A);
        double t2 = (-B + mysqrt(delta)) / (2 * A);
        if (dcmp(t1 - 1) <= 0 && dcmp(t1) >= 0)
            ret[num++] = point(x1 + t1 * dx, y1 + t1 * dy);
        if (dcmp(t2 - 1) <= 0 && dcmp(t2) >= 0)
            ret[num++] = point(x1 + t2 * dx, y1 + t2 * dy);
    }
}

// 取模
double abs(const point &o)
{
    return sqrt(dot(o, o));
}
// 定义叉积
double cross(const point &a, const point &b)
{
    return a.x * b.y - b.x * a.y;
}
// 求交点
point crosspt(const point &a, const point &b, const point &p, const point &q)
{
    double a1 = cross(b - a, p - a);
    double a2 = cross(b - a, q - a);
    return (p * a2 - q * a1) / (a2 - a1);
}

const int MAXN = 100010;
const double PI = acos(-1);
point res[MAXN];
double r;
// 圆与多边形交的面积
// 计算扇形的面积
double sector_area(const point &a, const point &b)
{
    double theta = atan2(a.y, a.x) - atan2(b.y, b.x);
    while (theta <= 0)
        theta += 2 * PI;
    while (theta > 2 * PI)
        theta -= 2 * PI;
    theta = min(theta, 2 * PI - theta);
    return r * r * theta / 2.0;
}
// 计算每一小三角形和圆面积交的面积
double calc(const point &a, const point &b)
{
    point p[2];
    int num = 0;
    int ina = dcmp(abs(a) - r) < 0;
    int inb = dcmp(abs(b) - r) < 0;
    if (ina)
    {
        if (inb)
            return fabs(cross(a, b)) / 2.0;
        else
        {
            circle_cross_segment(a, b, point(0, 0), r, p, num);
            return sector_area(b, p[0]) + fabs(cross(a, p[0])) / 2.0;
        }
    }
    else
    {
        if (inb)
        {
            circle_cross_segment(a, b, point(0, 0), r, p, num);
            return sector_area(p[0], a) + fabs(cross(p[0], b)) / 2.0;
        }
        else
        {
            circle_cross_segment(a, b, point(0, 0), r, p, num);
            if (num == 2)
                return sector_area(a, p[0]) + sector_area(p[1], b) + fabs(cross(p[0], p[1])) / 2.0;
            else
                return sector_area(a, b);
        }
    }
}
// 圆和多边形面积交
double area()
{
    double ret = 0;
    for (int i = 0; i < n; i++)
    {
        int sgn = dcmp(cross(res[i], res[(i + 1) % n]));
        if (sgn != 0)
            ret += sgn * calc(res[i], res[(i + 1) % n]);
    }
    return ret;
}

// 最小圆覆盖
// 三角形的外心
void circle_center(point p0, point p1, point p2, point &cp)
{
    point a = p1 - p0, b = p2 - p0;
    double c1 = (a.x * a.x + a.y * a.y) / 2, c2 = (b.x * b.x + b.y * b.y) / 2;
    double d = cross(a, b);
    cp.x = p0.x + (c1 * b.y - c2 * a.y) / d;
    cp.y = p0.y + (a.x * c2 - b.x * c1) / d;
}

// 两点求圆心(即中点坐标)
void circle_center(point p0, point p1, point &cp)
{
    cp = (p0 + p1) / 2;
}

point center;
double radius;
// 判断点是否在园内
bool point_in(const point &p)
{
    return cmp(dist(p, center) - radius) < 0;
}

// 最小圆覆盖
void min_circle_cover(point a[], int n)
{
    random_shuffle(a, a + n);
    radius = 0;
    center = a[0];
    for (int i = 1; i < n; i++)
        if (!point_in(a[i]))
        {
            center = a[i];
            radius = 0;
            for (int j = 0; j < i; j++)
                if (!point_in(a[j]))
                {
                    circle_center(a[i], a[j], center);
                    radius = dist(a[j], center);
                    for (int k = 0; k < j; k++)
                        if (!point_in(a[k]))
                        {
                            circle_center(a[i], a[j], a[k], center);
                            radius = dist(a[k], center);
                        }
                }
        }
}

point rotate(const point &p, double cost, double sint)
{
    double x = p.x, y = p.y;
    return point(x * cost - y * sint, x * sint + y * cost);
}

// 圆与圆求交
pair<point, point> crosspoint(point ap, double ar, point bp, double br)
{
    double d = (ap - bp).norm();
    double cost = (ar * ar + d * d - br * br) / (2 * ar * d);
    double sint = sqrt(1.0 - cost * cost);
    point v = (bp - ap) / (bp - ap).norm() * ar;
    return make_pair(ap + rotate(v, cost, -sint), ap + rotate(v, cost, sint));
}

// 圆的离散化
const int maxN = maxn * maxn + 3 * maxn;
struct Tcir
{
    double r;
    point o;
    void read()
    {
        scanf("%lf%lf%lf%lf", &o.x, &o.y, &r);
    }
};
struct Tinterval
{
    double x, y, Area, mid;
    int type;
    Tcir owner;
    void area(double l, double r)
    {
        double len = sqrt(sqr(l - r) + sqr(x - y));
        double d = sqrt(sqr(owner.r) - sqr(len) / 4.0);
        double angle = atan(len / 2.0 / d);
        Area = fabs(angle * sqr(owner.r) - d * len / 2.0);
    }
} inter[maxn];
double x[maxN], l;
int m, M, Mm;
bool compR(const Tcir &a, const Tcir &b)
{
    return a.r > b.r;
}
void Get(Tcir owner, double x, double &l, double &r)
{
    double y = fabs(owner.o.x - x);
    double d = sqrt(fabs(sqr(owner.r)) - sqr(y));
    l = owner.o.y + d, r = owner.o.y - d;
}
void Get_Interval(Tcir owner, double l, double r)
{
    Get(owner, l, inter[Mm].x, inter[Mm + 1].x);
    Get(owner, r, inter[Mm].y, inter[Mm + 1].y);
    Get(owner, (l + r) / 2.0, inter[Mm].mid, inter[Mm + 1].mid);
    inter[Mm].owner = inter[Mm + 1].owner = owner;
    inter[Mm].area(l, r);
    inter[Mm + 1].area(l, r);
    inter[Mm].type = 1;
    inter[Mm + 1].type = -1;
    Mm += 2;
}
bool comp_circle(const Tinterval &a, const Tinterval &b)
{
    return a.mid > b.mid + eps;
}
void Add(double xx)
{
    x[M++] = xx;
}
double dist2(const point &a, const point &b)
{
    return sqrt(dist(a, b));
}

// 圆的面积并1
double getUnion(int n, Tcir a[])
{
    int p = 0;
    sort(a, a + n, compR);
    for (int i = 0; i < m; i++)
    {
        bool f1 = true;
        for (int j = 0; j < i; j++)
            if (dist(a[i].o, a[j].o) <= sqr(a[i].r - a[j].r) + 1e-12)
            {
                f1 = false;
                break;
            }
        if (f1)
            a[p++] = a[i];
    }
    m = p;

    M = 0;
    for (int i = 0; i < m; i++)
    {
        Add(a[i].o.x - a[i].r);
        Add(a[i].o.x + a[i].r);
        Add(a[i].o.x);
        for (int j = i + 1; j < n; j++)
            if (dist2(a[i].o, a[j].o) <= sqr(a[i].r + a[j].r) + eps)
            {
                pair<point, point> cross = crosspoint(a[i].o, a[i].r, a[j].o, a[j].r);
                Add(cross.first.x);
                Add(cross.second.x);
            }
    }

    sort(x, x + M);
    p = 0;
    for (int i = 0; i < M; i++)
        if (!i || fabs(x[i] - x[i - 1]) > eps)
            x[p++] = x[i];
    M = p;

    double ans = 0;
    for (int i = 0; i < M; i++)
    {
        l = x[i], r = x[i + 1];
        Mm = 0;
        for (int j = 0; j < n; j++)
            if (fabs(a[j].o.x - l) < a[j].r + eps && fabs(a[j].o.x - r) < a[j].r + eps)
                Get_Interval(a[j], l, r);
        if (Mm)
        {
            sort(inter, inter + Mm, comp_circle);
            int cnt = 0;
            for (int i = 0; i < Mm; i++)
            {
                if (cnt > 0)
                {
                    ans += (fabs(inter[i - 1].x - inter[i].x) + fabs(inter[i - 1].y - inter[i].y)) * (r - l) / 2.0;
                    ans += inter[i - 1].type * inter[i - 1].Area;
                    ans -= inter[i].type * inter[i].Area;
                }
                cnt += inter[i].type;
            }
        }
    }
    return ans;
}

// 圆的面积并2
struct Circle
{
    point p;
    double r;

    bool operator<(const Circle &o) const
    {
        if (dcmp(r - o.r) != 0)
            return dcmp(r - o.r) == -1;
        if (dcmp(p.x - o.p.x) != 0)
            return dcmp(p.x - o.p.x) == -1;

        return dcmp(p.y - o.p.y) == 0;
    }

    bool operator==(const Circle &o) const
    {
        return dcmp(r - o.r) == 0 && dcmp(p.x - o.p.x) == 0 && dcmp(p.y - o.p.y) == 0;
    }
};
inline pair<point, point> crosspoint(const Circle &a, const Circle &b)
{
    return crosspoint(a.p, a.r, b.p, b.r);
}
Circle c[1000], tc[1000];
int cn, cm;
struct Node
{
    point p;
    double a;
    int d;

    Node(const point &p, double a, int d) : p(p), a(a), d(d) {}

    bool operator<(const Node &o) const
    {
        return a < o.a;
    }
};
double arg(point p)
{
    return arg(complex<double>(p.x, p.y));
}
double get_Union()
{
    sort(tc, tc + cm);
    cm = unique(tc, tc + cm) - tc;
    for (int i = cm - 1; i >= 0; i--)
    {
        bool ok = true;
        for (int j = i + 1; j < m; j++)
        {
            double d = (tc[i].p - tc[j].p).norm();
            if (dcmp(d - abs(tc[i].r - tc[j].r)) <= 0)
            {
                ok = false;
                break;
            }
        }
        if (ok)
            c[cn++] = tc[i];
    }
    double ans = 0;
    for (int i = 0; i < n; i++)
    {
        vector<Node> event;
        point boundary = c[i].p + point(-c[i].r, 0);
        event.push_back(Node(boundary, -PI, 0));
        event.push_back(Node(boundary, PI, 0));
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            double d = (c[i].p - c[j].p).norm();
            if (dcmp(d - (c[i].r + c[j].r)) < 0)
            {
                pair<point, point> ret = crosspoint(c[i], c[j]);
                double x = arg(ret.first - c[i].p);
                double y = arg(ret.second - c[i].p);
                if (dcmp(x - y) > 0)
                {
                    event.push_back(Node(ret.first, x, 1));
                    event.push_back(Node(boundary, PI, -1));
                    event.push_back(Node(boundary, -PI, 1));
                    event.push_back(Node(ret.second, y, -1));
                }
                else
                {
                    event.push_back(Node(ret.first, x, 1));
                    event.push_back(Node(ret.second, y, -1));
                }
            }
        }

        sort(event.begin(), event.end());
        int sum = event[0].d;
        for (int j = 1; j < (int)event.size(); j++)
        {
            if (sum == 0)
            {
                ans += cross(event[j - 1].p, event[j].p);
                double x = event[j - 1].a;
                double y = event[j].a;
                double area = c[i].r * c[i].r * (y - x) / 2.;
                point v1 = event[j - 1].p - c[i].p;
                point v2 = event[j].p - c[i].p;
                area -= cross(v1, v2) / 2.;
                ans += area;
            }
            sum += event[j].d;
        }
    }
    return ans;
}

void solve()
{
}
int main()
{
    // ios::sync_with_stdio(false);
    // cin.tie(0), cout.tie(0);
    // cout << fixed << setprecision(10) << area() << '\n';
    int T;
    cin >> T;
    while (T--)
        solve();
    system("pause");
    return 0;
}