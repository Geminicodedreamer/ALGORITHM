// 4177. 曲线
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;

const int N = 100010;
const double eps = 1e-9;

int a[N], b[N], c[N], n;
double check(double x)
{
    double res = -1e9;

    for (int i = 0; i < n; i++)
        res = max(res, a[i] * x * x + b[i] * x + c[i]);

    return res;
}
int main()
{
    int T;
    scanf("%d", &T);
    while (T--)
    {
        scanf("%d", &n);
        for (int i = 0; i < n; i++)
            scanf("%d%d%d", &a[i], &b[i], &c[i]);
        double l = 0.0, r = 1000.0;
        while (r - l > eps)
        {
            double d = (r - l) / 3.0;
            double lmid = l + d, rmid = r - d;
            if (check(lmid) > check(rmid))
                l = lmid;
            else
                r = rmid;
        }
        printf("%.4lf\n", check(l));
    }
}