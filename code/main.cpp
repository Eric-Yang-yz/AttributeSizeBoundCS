#include"utility.h"
#include"Snap.h"
#include"main.h"

int main(int argc,char *argv[])
{
	PUNGraph G_;
	TStr path = "C:\\¿ÆÑÐ\\code_yyz\\data_set\\test.graph.yyz.txt";
	printf("Start Loading Network...\n");
	G_ = LoadNetwork(path);
	printf("Finish Loading Network...\n");
	G = G_;

	Solve(2, 7, 11);
	if (k_hat == 0) {
		printf("feasible H not found!\n");
		return 0;
	}
	printf("size: %d, k: %d\n", H->GetNodes(), k_hat);
	printf("vertices in H:\n");
	for (auto it = H->BegNI(); it < H->EndNI(); it++) {
		printf("%d\n", it.GetId());
	}
	printf("edges in H:\n");
	for (auto it = H->BegEI(); it < H->EndEI(); it++) {
		printf("%d %d\n", it.GetSrcNId(), it.GetDstNId());
	}
}