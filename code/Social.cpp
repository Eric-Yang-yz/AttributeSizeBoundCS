#include"Social.h"
using namespace std;

PUNGraph LoadNetwork(TStr path)
{
	auto G = TSnap::LoadEdgeList<PUNGraph>(path, 0, 1);
	return G;
}