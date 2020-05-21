#include<iostream>
#include <fstream>
#include<ctime>
#include<random>
/*

 !!! В работе был использован Муравьиный (ACO TSP) алгоритм, описанный в видео https://www.youtube.com/watch?v=EwDP_bAb-OI&feature=emb_logo 
 !!! Вся терминалогия так же взята от туда

*/


using namespace std;
int N;
double tau[1000][1000];
double Lij[1000][1000];
bool badWay[1000];
double probability[2][1000];
int way[1000][1000];
double wayLen[1000];

uniform_real_distribution<double> getRandom(0, 100);
mt19937 randGen(time(NULL));

void wayACO(int i, int taken, bool topAnt, int noVisited);
void ACO();

void ACO() {
	int minL = 11000000;
	for (int iter = 0; iter < 200; iter++) {
		bool topAnt = false;
		for (int taken = 1; taken <= N; taken++) {
			badWay[taken] = true;
			way[taken][0] = taken;
			wayACO(taken, taken, topAnt, N - 1);

			if (minL > wayLen[taken]) minL = wayLen[taken];
			for (int j = 1; j <= N; j++) badWay[j] = 0;
		}

		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++) tau[i][j] /= 2;

		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++) tau[way[i][j - 1]][way[i][j]] += (1.0 / wayLen[i]);

		for (int i = 0; i <= N; i++) {
			for (int j = 0; j <= N; j++) way[i][j] = 0;
			wayLen[i] = 0;
		}
	}

	cout << minL;
}

void wayACO(int i, int taken, bool topAnt, int noVisited) {
	double sum = 0.0;
	for (int j = 1; j <= N; j++)
		if (!badWay[j]) sum += pow(1.0 / Lij[i][j], 2.0) * tau[i][j];

	int probN = 0;
	for (int j = 1; j <= N; j++) {
		if (!badWay[j]) {
			probability[0][probN] = 100 * pow(1.0 / Lij[i][j], 2.0) * tau[i][j] / sum;
			probability[1][probN] = j;
			probN++;
		}
	}

	if (probN == 1) {
		way[taken][N - noVisited] = probability[1][0];
		int w = (int)probability[1][0];
		wayLen[taken] += Lij[i][w];
		badWay[w] = true;
		wayLen[taken] += Lij[w][taken];
		way[taken][N] = taken;
		return;
	}

	if (topAnt) {
		double mostProb = 0.0;
		int k;
		for (int j = 0; j < probN; j++) {
			if (probability[0][j] > mostProb) {
				mostProb = probability[0][j];
				k = probability[1][j];
			}
		}
		wayLen[taken] += Lij[i][k];
		badWay[k] = true;
		wayACO(k, taken, topAnt, noVisited - 1);
		return;
	}
	else {
		for (int j = 1; j < probN; j++) {
			probability[0][j] += probability[0][j - 1];
		}

		double choice = getRandom(randGen);

		for (int j = 0; j < probN; j++) {
			if (choice < probability[0][j]) {
				way[taken][N - noVisited] = probability[1][j];
				int w = (int)probability[1][j];
				wayLen[taken] += Lij[i][w];
				badWay[w] = true;
				wayACO(w, taken, topAnt, noVisited - 1);
				return;
			}
		}
	}
}


int main(void) {
	ifstream cin;
	//cin.open("test/five_d.txt");
	//cin.open("test/p01_d.txt");
	//cin.open("test/gr17_d.txt");
	//cin.open("test/fri26_d.txt");
	cin.open("test/att48_d.txt");

	cin >> N;
	for (int i = 1; i <= N; i++)
		for (int j = 1; j <= N; j++) {
			cin >> Lij[i][j];
			tau[i][j] = getRandom(randGen);
		}

	ACO();

	return 0;
}
