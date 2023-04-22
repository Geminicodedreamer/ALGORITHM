// 4799. 最远距离
#include <iostream>
#include <cstring>
#include <algorithm>

using namespace std;
const int N = 100010, M = 2 * N;
int e[M], ne[M], h[N], idx;
typedef pair<int, int> PII;

void add(int a, int b) // 添加一条边a->b
{
    e[idx] = b, ne[idx] = h[a], h[a] = idx++;
}

PII dfs(int u, int st)
{
    PII res = {0, u};
    for (int i = h[u]; ~i; i = ne[i])
    {
        auto ver = e[i];
        if (ver == st)
            continue;
        auto t = dfs(ver, u);
        t.first++;
        if (res < t)
            res = t;
    }
    return res;
}

int main()
{
    int n, m;
    scanf("%d%d", &n, &m);
    memset(h, -1, sizeof h);
    while (m--)
    {
        int a, b;
        scanf("%d%d", &a, &b);
        add(a, b), add(b, a);
    }

    int a = dfs(1, -1).second;
    cout << dfs(a, -1).first << '\n';
}
