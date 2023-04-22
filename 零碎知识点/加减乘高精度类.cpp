#include<bits/stdc++.h>
#define ull unsigned long long
using namespace std;
const ull L = 23333, A = 9, N = 1e9;
struct super_int {
	ull l;
	ull d[L] = {0};
};
super_int operator + (const super_int& a, const super_int& b) {
	super_int c;
	for (ull i = 0, s = 0; i < L; i++) {
		s = (c.d[i] = s + a.d[i] + b.d[i]) / N;
		c.d[i] %= N;
	}
	return c;
}
super_int operator - (const super_int& a) {
	super_int si1, b;
	si1.d[0] = 1;
	for (int i = 0; i < L; i++)
		b.d[i] = N - 1 - a.d[i];
	return b + si1;
}
super_int operator - (const super_int& a, const super_int& b) {
	return a + -b;
}
super_int operator * (const super_int& a, const super_int& b) {
	super_int c;
	for (ull k = 0, s = 0; k < L; k++) {
		c.d[k] = s;
		s = 0;
		for (int i = 0; i <= k; i++) {
			int j = k - i;
			s += (c.d[k] += a.d[i] * b.d[j]) / N;
			c.d[k] %= N;
		}
	}
	return c;
}
istream& operator >> (istream& in, super_int& a) {
	string s;
	in >> s;
	bool sgn = s[0] == '-';
	for (ull i = sgn, di = 0; i < s.length(); i++) {
		di = di * 10 + s[i] - '0';
		if ( (s.length() - i) % A == 1) {
			a.d[ (s.length() - i) / A] = di;
			di = 0;
		}
	}
	if (sgn)  a = -a;
	return in;
}
ostream& operator << (ostream& out, const super_int& a) {
	super_int b = a;
	bool sgn = a.d[L - 1] == N - 1, flag = 0;
	if (sgn) {
		out << '-';
		b = -a;
	}
	for (int i = L - 1 - sgn; i >= 0; i--)
		if (flag)  out << setw (A) << setfill ('0') << b.d[i];
		else  if (b.d[i] || !i) {
			flag = 1;
			out << b.d[i];
		}
	return out;
}
int main() {
	super_int a, b;
	cin >> a >> b;
	cout << a + b;
	return 0;
}