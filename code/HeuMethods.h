#pragma once
#ifndef HEU_METHODS_H
#define HEU_METHODS_H
#include"SC-BRB.h"


void DegreeBasedReduction(TIntV &C, TIntV &R);
void DistanceBasedReduction(TIntV &C, TIntV &R);
int BinarySearchD(int k, int h);
int DegreeBasedUB(TIntV C, TIntV R);
int NeighborReconstructUB(TIntV C, TIntV R);
int OneQuickSort(vector<int>& a, int left, int right, int mode);
void FindLargestK(vector<int>& a, int k);
void FindSmallestK(vector<int>& a, int k);
int DegreeClassificationUB(TIntV C, TIntV R);
double GetConnectionScore(TIntV C, int v);
vector<int> GetDominateSet(TIntV C, TIntV R, int v);
#endif