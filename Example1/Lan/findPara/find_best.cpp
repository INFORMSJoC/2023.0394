#include <bits/stdc++.h>
using namespace std;
double width[10];
string num2str(int x) {
	stringstream ss;
	ss << x;
	return ss.str();
}
int main() {
	int replications = 20;
	for (int i = 1; i < 10; ++ i) {
		ifstream fin;
		string directory_name = "group" + num2str(i);
		directory_name = directory_name + "/results/";
		width[i] = 0.0;
		for (int k = 1; k <= 20; ++ k) {
			string filename = directory_name + num2str(k) + ".txt";
			fin.open(filename);
			double ll, rr;
			fin >> ll >> rr;
			fin.close();
			width[i] += rr - ll;
		}
		width[i] /= 1.0 * replications;
	}
	ofstream fout;
	fout.open("summary.txt");
	double min_v = 1e6;
	int min_id = -1;
	for (int k = 1; k < 10; ++ k) fout << width[k] << endl;
	for (int k = 1; k < 10; ++ k)
		if (min_v > width[k]) {
			min_v = width[k];
			min_id = k;
		}
	fout << min_id << endl;
	fout.close();
	return 0;
}
