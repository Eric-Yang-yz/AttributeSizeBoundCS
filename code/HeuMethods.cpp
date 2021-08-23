#include"HeuMethods.h"
const int INF = 0x7fffffff;

void DegreeBasedReduction(TIntV &C, TIntV &R)
{
	int d_CR;		//refer to d(C U R) v
	int d_Cv;		//refer to d(C U v) v + High - |C| - 1;
	TIntV RemoveSet;

	for (auto it = R.BegI(); it < R.EndI(); it++) {
		int v = it->Val;
		PUNGraph subG;
		TIntV SubVertice = C;
		/* Calculate d_Cv */
		SubVertice.Add(v);
		subG = TSnap::GetSubGraph(G, SubVertice);
		d_Cv = subG->GetNI(v).GetDeg() + high - C.Len() - 1;
		/* Calculate d_CR */
		SubVertice.DelLast();
		SubVertice.AddV(R);
		subG = TSnap::GetSubGraph(G, SubVertice);
		d_CR = subG->GetNI(v).GetDeg();

		subG.Clr();
		if (min(d_CR, d_Cv) <= k_hat)
			RemoveSet.Add(v);
	}
	for (auto it = RemoveSet.BegI(); it < RemoveSet.EndI(); it++) {
		R.DelIfIn(it->Val);
	}
}

void DistanceBasedReduction(TIntV& C, TIntV& R)
{
	int maxD = BinarySearchD(k_hat + 1, high);
	TIntV SubVertice = C;
	PUNGraph subG;
	TIntV RemoveSet;
	for (auto it = R.BegI(); it < R.EndI(); it++) {
		SubVertice.Add(it->Val);
		subG = TSnap::GetSubGraph(G, SubVertice);

		/* start bfs to find max distance */
		int dis = 0;
		unordered_map<int, bool> visit;
		for (int i = 0; i < SubVertice.Len(); i++) {
			visit[i] = 0;
		}
		queue<pair<int,int> > qu;			//pair<nodeID, distance>
		qu.push(make_pair(it->Val,0));
		while (!qu.empty()) {
			pair<int, int> tmp = qu.front();
			qu.pop();
			if (visit[tmp.first])continue;
			visit[tmp.first] = 1;
			dis = max(dis, tmp.second);
			auto sub_it = subG->GetNI(tmp.first);
			for (int i = 0; i < sub_it.GetDeg(); i++) {
				int nbr = sub_it.GetNbrNId(i);
				if (!visit[nbr]) {
					qu.push(make_pair(nbr, tmp.second + 1));
				}
			}
		}
		if (dis > maxD)RemoveSet.Add(it->Val);
		
		subG.Clr();
		SubVertice.DelLast();
	}
	for (auto it = RemoveSet.BegI(); it < RemoveSet.EndI(); it++) {
		R.DelIfIn(it->Val);
	}
}

int BinarySearchD(int k,int h)
{
	if (k == 1)return high - 1;
	int left, right, maxD;
	left = 1; right = G->GetNodes();
	while (left <= right) {
		int val;
		int mid = left + right >> 1;
		if (1 <= mid && mid <= 2)val = k + mid;
		else val = k + mid + 1 + (mid / 3) * (k - 2);
		if (val > high)right = mid - 1;
		else if (val == high) {
			maxD = mid;
			break;
		}
		else {
			maxD = mid;
			left = mid + 1;
		}
	}
	return maxD;
}

int DegreeBasedUB(TIntV C, TIntV R)
{
	int ub = G->GetNI(TSnap::GetMxDegNId(G)).GetDeg();
	TIntV CR = C;
	CR.AddV(R);
	PUNGraph subG = TSnap::GetSubGraph(G, CR);
	PUNGraph commC = TSnap::GetSubGraph(G, C);
	for (auto it = C.BegI(); it < C.EndI(); it++) {
		int d_CR;
		int d_Cu;
		d_CR = subG->GetNI(it->Val).GetDeg();
		d_Cu = commC->GetNI(it->Val).GetDeg() + high - C.Len();
		ub = min(ub, min(d_CR, d_Cu));
	}
	subG.Clr();
	commC.Clr();
	return ub;
}

int NeighborReconstructUB(TIntV C, TIntV R)
{
	int ub;
	TIntV SubVertice;
	PUNGraph subG;
	vector<int> Rset;		//refer to vertices in R adjcent to C
	vector<int> Rdegree;	//refer to each degree in R'  (CU{v})
	vector<int> Cdegree;	//refer to each degree in C
	SubVertice = R;
	/* find vertices in R which are adjacent to C */
	for (auto it = C.BegI(); it < C.EndI(); it++) {
		SubVertice.Add(it->Val);
		subG = TSnap::GetSubGraph(G, SubVertice);
		auto curNode = subG->GetNI(it->Val);
		for (int i = 0; i < curNode.GetDeg(); i++) {
			Rset.push_back(curNode.GetNbrNId(i));
		}
		SubVertice.DelLast();
		subG.Clr();
	}
	// Remove duplicate
	auto ed = unique(Rset.begin(), Rset.end());
	Rset.erase(ed, Rset.end());
	/* find R' */
	SubVertice = C;
	for (auto ii : Rset) {
		SubVertice.Add(ii);
		subG = TSnap::GetSubGraph(G, SubVertice);
		Rdegree.push_back(subG->GetNI(ii).GetDeg());
		SubVertice.DelLast();
		subG.Clr();
	}
	FindLargestK(Rdegree, high - C.Len());
	/* get Cdegree */
	subG = TSnap::GetSubGraph(G, SubVertice);
	for (auto it = C.BegI(); it < C.EndI(); it++) {
		Cdegree.push_back(subG->GetNI(it->Val).GetDeg());
	}
	subG.Clr();
	/* calculate upper bound */
	for (int i = 0; i < min(Rdegree.size(),high - C.Len()); i++) {
		FindSmallestK(Cdegree, Rdegree[i]);
		for (int j = 0; j < Rdegree[i]; j++)Cdegree[j]++;
	}
	ub = Cdegree[0];
	for (int i = 1; i < Cdegree.size(); i++)
		ub = min(ub, Cdegree[i]);
	return ub;
}

/* implement quick sort one turn,returns the index of the middle element */
int OneQuickSort(vector<int>& a, int left, int right,int mode)
{
	int i, j, tmp;
	i = left; j = right;
	tmp = a[left];
	if (mode) {		//large->small
		while (1) {
			while (i < j && a[j] < tmp)j--;
			if (i == j) {
				a[i] = tmp;
				break;
			}
			else {
				a[i] = a[j];
				i++;
			}
			while (i<j && a[i]>tmp)i++;
			if (i == j) {
				a[i] = tmp;
				break;
			}
			else {
				a[j] = a[i];
				j--;
			}
		}
	}
	else {			//small->large
		while (1) {
			while (i < j && a[j] > tmp)j--;
			if (i == j) {
				a[i] = tmp;
				break;
			}
			else {
				a[i] = a[j];
				i++;
			}
			while (i<j && a[i]<tmp)i++;
			if (i == j) {
				a[i] = tmp;
				break;
			}
			else {
				a[j] = a[i];
				j--;
			}
		}
	}
	return i;
}
/* find largest k elements in a list with complexity O(n) */
void FindLargestK(vector<int> &a, int k)
{
	if (k >= a.size() || k == 0)return;
	int left, right;
	left = 0; right = a.size() - 1;
	while (left < right) {
		int idx = OneQuickSort(a, left, right, 1);
		//index start from 0
		//if the middle element has index k-1 or k, then we can ensure the top k elements
		if (idx-left == k - 1 || idx-left == k) {
			break;
		}
		else if (idx-left > k) {
			right = idx - 1;
		}
		else {
			k -= idx - left + 1;		//update k
			left = idx + 1;
		}
	}
	
}
/* find smallest k elements in a list with complexity O(n) */
void FindSmallestK(vector<int>& a, int k)
{
	if (k >= a.size() || k == 0)return;
	int left, right;
	left = 0; right = a.size() - 1;
	while (left < right) {
		int idx = OneQuickSort(a, left, right, 0);
		//index start from 0
		//if the middle element has index k-1 or k, then we can ensure the top k elements
		if (idx - left == k - 1 || idx - left == k) {
			break;
		}
		else if (idx - left > k) {
			right = idx - 1;
		}
		else {
			k -= idx - left + 1;		//update k
			left = idx + 1;
		}
	}
}

int DegreeClassificationUB(TIntV C, TIntV R)
{
	int ub = INF;
	TIntV SubVertice;
	PUNGraph subG;
	vector<pair<int,int> > Cdegree;		//pair<nodeID,degree> in subG induced by C
	TIntV Ctao;

	int dminC = INF, dmaxC = -1;
	subG = TSnap::GetSubGraph(G, C);
	for (auto it = subG->BegNI(); it < subG->EndNI(); it++) {
		dminC = min(dminC, it.GetDeg());
		dmaxC = max(dmaxC, it.GetDeg());
		Cdegree.push_back(make_pair(it.GetId(), it.GetDeg()));
	}
	subG.Clr();
	for (int tao = dminC; tao <= dmaxC; tao++) {
		for (auto ii : Cdegree) {
			if (ii.second == tao)Ctao.Add(ii.first);
		}
		/* find vertices in R which are adjacent to Ctao */
		vector<int> Rset;
		vector<int> Rdegree;
		vector<int> Ctaodegree;
		SubVertice = R;
		for (auto it = Ctao.BegI(); it < Ctao.EndI(); it++) {
			SubVertice.Add(it->Val);
			subG = TSnap::GetSubGraph(G, SubVertice);
			auto curNode = subG->GetNI(it->Val);
			for (int i = 0; i < curNode.GetDeg(); i++) {
				Rset.push_back(curNode.GetNbrNId(i));
			}
			SubVertice.DelLast();
			subG.Clr();
		}
		//Remove duplicate
		auto ed = unique(Rset.begin(), Rset.end());
		Rset.erase(ed, Rset.end());
		/* find R' */
		SubVertice = Ctao;
		for (auto ii : Rset) {
			SubVertice.Add(ii);
			subG = TSnap::GetSubGraph(G, SubVertice);
			Rdegree.push_back(subG->GetNI(ii).GetDeg());
			SubVertice.DelLast();
			subG.Clr();
		}
		FindLargestK(Rdegree, high - C.Len());
		/* get Ctaodegree */
		subG = TSnap::GetSubGraph(G, SubVertice);
		for (auto it = Ctao.BegI(); it < Ctao.EndI(); it++) {
			Ctaodegree.push_back(subG->GetNI(it->Val).GetDeg());
		}
		subG.Clr();
		/* calculate upper bound */
		for (int i = 0; i < min(Rdegree.size(), high - C.Len()); i++) {
			FindSmallestK(Ctaodegree, Rdegree[i]);
			for (int j = 0; j < Rdegree[i]; j++)Ctaodegree[j]++;
		}
		int ubtao = Ctaodegree[0];
		for (int i = 1; i < Ctaodegree.size(); i++) {
			ubtao = min(ubtao, Ctaodegree[i]);
		}
		ub = min(ub, ubtao);
		if (ub <= tao + 1)break;
	}
	return ub;
}

double GetConnectionScore(TIntV C, int v)
{
	double score = 0;
	PUNGraph subG,subGc;
	TIntV SubVertice;
	SubVertice = C;
	SubVertice.Add(v);
	subGc = TSnap::GetSubGraph(G, C);
	subG = TSnap::GetSubGraph(G, SubVertice);
	auto curNode = subG->GetNI(v);
	for (int i = 0; i < curNode.GetDeg(); i++) {
		auto c_node = subGc->GetNI(curNode.GetNbrNId(i));
		score += 1.0 / (c_node.GetDeg() * 1.0);
	}
	subG.Clr(); subGc.Clr();
	return score;
}

vector<int> GetDominateSet(TIntV C, TIntV R, int v)
{
	vector<int> DominSet;
	TIntV SubVertice;
	PUNGraph subG;
	unordered_map<int, bool> isR, isNbrV;

	for (auto it = C.BegI(); it < C.EndI(); it++) {
		isR[it->Val] = 0;
		isNbrV[it->Val] = 0;
	}
	for (auto it = R.BegI(); it < R.EndI(); it++) {
		isR[it->Val] = 1;
		isNbrV[it->Val] = 0;
	}
	SubVertice = C;
	SubVertice.AddV(R);
	subG = TSnap::GetSubGraph(G, SubVertice);
	auto nodeV = subG->GetNI(v);
	isNbrV[v] = 1;
	for (int i = 0; i < nodeV.GetDeg(); i++) {
		isNbrV[nodeV.GetNbrNId(i)] = 1;
	}
	for (int i = 0; i < nodeV.GetDeg(); i++) {
		int cur = nodeV.GetNbrNId(i);
		if (!isR[cur])continue;
		auto it_cur = subG->GetNI(cur);
		bool flag = true;
		for (int j = 0; j < it_cur.GetDeg(); j++) {
			if (!isNbrV[it_cur.GetNbrNId(j)]) {
				flag = false;
				break;
			}
		}
		if (flag)DominSet.push_back(cur);
	}
	subG.Clr();
	return DominSet;
}