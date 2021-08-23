#include"SC-BRB.h"
PUNGraph G;		//total graph
PUNGraph H;		//the current best community
int k_hat;		//current biggest k, 0 means not found yet
int low, high;

bool cmp(pair<int, double> a, pair<int, double> b)
{
	return a.second > b.second;
}

void BRB(TIntV C, TIntV R)
{
	DegreeBasedReduction(C, R);
	DistanceBasedReduction(C, R);
	int k_cur = GetMinDegree(C);
	if (low <= C.Len() && C.Len() <= high && k_cur > k_hat) {
		k_hat = k_cur;
		H = TSnap::GetSubGraph(G, C);
	}
	if (DegreeBasedUB(C, R) <= k_hat)return;
	if (NeighborReconstructUB(C, R) <= k_hat)return;
	if (DegreeClassificationUB(C, R) <= k_hat)return;
	if (C.Len() < high && R.Len() > 0) {
		int Vstar = -1;
		double maxscore = 0;
		unordered_map<int, double> allscore;				//<nodeID,ConnectionScore>
		for (auto it = R.BegI(); it < R.EndI(); it++) {
			double tscore = GetConnectionScore(C, it->Val);
			allscore[it->Val] = tscore;
			if (tscore > maxscore) {
				maxscore = tscore;
				Vstar = it->Val;
			}
		}
		if (Vstar != -1) {
			vector<int> dominSet = GetDominateSet(C, R, Vstar);
			vector<pair<int, double> > dominSet_s;
			for (auto ii : dominSet) {
				dominSet_s.push_back(make_pair(ii, allscore[ii]));
			}
			sort(dominSet_s.begin(), dominSet_s.end(), cmp);
			for (int i = 0; i < dominSet_s.size(); i++) {
				TIntV Cnew = C;
				TIntV Rnew = R;
				Cnew.Add(Vstar); Cnew.Add(dominSet_s[i].first);
				Rnew.DelIfIn(Vstar);
				for (int j = 0; j <= i; j++) {
					Rnew.DelIfIn(dominSet_s[j].first);
				}
				BRB(Cnew, Rnew);
			}
			TIntV Cnew = C;
			TIntV Rnew = R;
			Cnew.Add(Vstar); Rnew.DelIfIn(Vstar);
			for (int i = 0; i < dominSet.size(); i++) {
				Rnew.DelIfIn(dominSet[i]);
			}
			BRB(Cnew, Rnew);
			Cnew.DelLast();
			BRB(Cnew, Rnew);
		}
	}
}
/* initialize variables and then start BRB */
void Solve(int query,int l,int h)
{
	/* Init */
	TIntV C, R;
	C.Add(query);
	for (auto it = G->BegNI(); it < G->EndNI(); it++) {
		if (it.GetId() == query)continue;
		else R.Add(it.GetId());
	}
	low = l;
	high = h;
	k_hat = 0;
	
	/* start BRB */
	BRB(C, R);
}

int GetMinDegree(TIntV C)
{
	int mindeg = 0x7fffffff;
	PUNGraph subG = TSnap::GetSubGraph(G, C);
	for (auto it = subG->BegNI(); it < subG->EndNI(); it++) {
		mindeg = min(mindeg, it.GetDeg());
	}
	return mindeg;
}