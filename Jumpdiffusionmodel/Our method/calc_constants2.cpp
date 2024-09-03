#include <bits/stdc++.h>
using namespace std;
const int n = 10;
const double m = 0.03, theta = 0.15, lambda = 0.25 / 50;
const double pi = acos(-1.0);
double getNp(double v, double mu, double sigma) {
	double _v = (v - mu) / sigma;
	return erfc(-_v / sqrt(2.0)) / 2.0;
}
double getNp(double ll, double rr, double mu, double sigma) {
	return getNp(rr, mu, sigma) - getNp(ll, mu, sigma);
}
vector<double> get_interval() {
	vector<double> ret; ret.clear();
	double last = -1e10;
	for (int i = 1; i < n; ++ i) {
		double v = 1.0 * i / n;
		double ll = last, rr = 1e10;
		while (rr - ll > 1e-8) {
			double mid = (rr + ll) / 2;
			if (getNp(mid, m, theta) >= v) rr = mid; else ll = mid;
		}
		ret.push_back(rr);
		last = rr;
		cout << rr << endl;
	}
	return ret;
}
void test_intervals(vector<double> intervals) {
	default_random_engine generator;
	normal_distribution<double> N01(0.0, 1.0);
	vector<int> cnt; cnt.clear();
	for (int i = 0; i < n; ++ i) cnt.push_back(0);
	int nn = 100000000;
	for (int i = 0; i < nn; ++ i) {
		double v = N01(generator) * theta + m;
		if (v < intervals[0]) cnt[0] = cnt[0] + 1; else {
			for (int j = n - 2; j >= 0; -- j)
				if (v >= intervals[j]) {
					cnt[j + 1] = cnt[j + 1] + 1;
					break;
				}
		}
	}
	for (int i = 0; i < n; ++ i) cout << 1.0 * cnt[i] / nn << ' '; cout << endl;
}
double calc_m1_std(double d, double u) {
	double sign = 1.0;
	if (u + d < 0.0) sign = -1.0;
	double uu = max(abs(u), abs(d)), vv = min(abs(u), abs(d));
	uu = uu * uu / 2; vv = vv * vv / 2;
	return (exp(-vv) - exp(-uu)) / sqrt(pi * 2) * sign;
}
double calc_m1(double ll, double rr, double mu, double sigma) {
	double u = (rr - mu) / sigma, d = (ll - mu) / sigma;
	return mu * getNp(d, u, 0.0, 1.0) + sigma * calc_m1_std(d, u);
}

double calc_m1(double ll, double rr) {
	return calc_m1(ll, rr, m, theta) * (1.0 - exp(-lambda));
}
vector<double> calc_m1(vector<double> intervals) {
	vector<double> ret; ret.clear();
	ret.push_back(calc_m1(-1e10, intervals[0]));
	for (int i = 0; i + 1 < int(intervals.size()); ++ i) 
		ret.push_back(calc_m1(intervals[i], intervals[i + 1]));
	ret.push_back(calc_m1(intervals[int(intervals.size()) - 1], 1e10));
	for (int i = 0; i < ret.size(); ++ i) cout << ret[i] << ' '; cout << endl;
	return ret;
}
void test_m1(vector<double> intervals) {
	default_random_engine generator;
	bernoulli_distribution distribution(1.0 - exp(-lambda));
	normal_distribution<double> N01(0.0, 1.0);
	vector<double> cnt; cnt.clear();
	for (int i = 0; i < n - 1; ++ i) cnt.push_back(0.0);
	int nn = 1000000000;
	for (int i = 0; i < nn; ++ i) {
		int t1 = distribution(generator);
		if (t1 == 0) continue;
		double mu = m * t1, sigma = theta * sqrt(1.0 * t1);
		double v = N01(generator) * theta + m;
		if (v < intervals[0]) cnt[0] = cnt[0] + v; else {
			for (int j = n - 2; j >= 0; -- j)
				if (v >= intervals[j]) {
					cnt[j + 1] = cnt[j + 1] + v;
					break;
				}
		}
	}
	for (int i = 0; i < n; ++ i) cout << cnt[i] / nn << ' '; cout << endl;
}

double calc_m2_std(double v) {
	double ret = -v * exp(-v * v / 2);
	ret += sqrt(pi / 2) * erf(v / sqrt(2.0));
	return ret / sqrt(pi * 2);
}
double calc_m2(double ll, double rr, double mu, double sigma) {
	double u = (rr - mu) / sigma, d = (ll - mu) / sigma;
	double uu = max(abs(u), abs(d)), vv = min(abs(u), abs(d));
	double ret = mu * mu * getNp(d, u, 0.0, 1.0) + calc_m1_std(d, u) * mu * sigma * 2;
	double tt;
	if (u * d <= 0.0)
		tt = calc_m2_std(abs(d)) + calc_m2_std(abs(u));
	else 
		tt = calc_m2_std(uu) - calc_m2_std(vv);
	ret += tt * sigma * sigma;
	return ret;
}
double calc_m2(double ll, double rr) {
	return calc_m2(ll, rr, m, theta) * (1.0 - exp(-lambda));
}
vector<double> calc_m2(vector<double> intervals) {
	vector<double> ret; ret.clear();
	ret.push_back(calc_m2(-1e10, intervals[0]));
	for (int i = 0; i + 1 < int(intervals.size()); ++ i) 
		ret.push_back(calc_m2(intervals[i], intervals[i + 1]));
	ret.push_back(calc_m2(intervals[int(intervals.size()) - 1], 1e10));
	for (int i = 0; i < ret.size(); ++ i) cout << ret[i] << ' '; cout << endl;
	return ret;
}
void test_m2(vector<double> intervals) {
	default_random_engine generator;
	bernoulli_distribution distribution(1.0 - exp(-lambda));
	normal_distribution<double> N01(0.0, 1.0);
	vector<double> cnt; cnt.clear();
	for (int i = 0; i < n - 1; ++ i) cnt.push_back(0.0);
	int nn = 1000000000;
	for (int i = 0; i < nn; ++ i) {
		int t1 = distribution(generator);
		if (t1 == 0) continue;
		double v = N01(generator) * theta + m;
		if (v <= intervals[0]) cnt[0] = cnt[0] + v * v; else {
			for (int j = n - 2; j >= 0; -- j)
				if (v > intervals[j]) {
					cnt[j + 1] = cnt[j + 1] + v * v;
					break;
				}
		}
	}
	for (int i = 0; i < n; ++ i) cout << cnt[i] / nn << ' '; cout << endl;
}

int main() {
	vector<double> intervals = get_interval();
	//test_intervals(intervals);
	vector<double> m1 = calc_m1(intervals);
	//test_m1(intervals);
	
	vector<double> m2 =	calc_m2(intervals);
	//test_m2(intervals);
	
	freopen("interval2.txt", "w", stdout);
	printf("left_points\n");
	printf("-1e10"); 
	for (int i = 0; i < intervals.size(); ++ i) printf(" %.8lf", intervals[i]);
	printf("\n");
	printf("right_points\n");
	for (int i = 0; i < intervals.size(); ++ i) printf("%.8lf ", intervals[i]);
	printf("1e10\n");
	printf("yk\n");
	
	for (int i = 0; i <= intervals.size(); ++ i)
	if (i < (n + 1) / 2) printf(" %.8lf", intervals[i]); else printf(" %.8lf", intervals[i - 1]);
	printf("\n");
	/*
	printf("%.8lf ", intervals[0] - (intervals[1] - intervals[0]) / 2); 
	for (int i = 0; i < intervals.size() - 1; ++ i)
		printf("%.8lf ", (intervals[i] + intervals[i + 1]) / 2);
	printf("%.8lf\n", intervals[intervals.size() - 1] + (intervals[intervals.size() - 1] - intervals[intervals.size() - 2]) / 2);
	*/
	printf("mean\n");
	for (int i = 0; i < m1.size(); ++ i) printf(" %.8lf", m1[i]);
	printf("\n");
	printf("variance\n");
	for (int i = 0; i < m1.size(); ++ i) printf(" %.8lf", m2[i] - m1[i] * m1[i]);
	printf("\n");
	printf("constant\n");
	for (int i = 0; i < m1.size(); ++ i) {
		double v = m2[i] - m1[i] * m1[i];
		v = v / m1[i];
		printf(" %.8lf", v);
	}
	printf("\n");
	
	return 0;
}
