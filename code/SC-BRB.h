#pragma once
#ifndef SCBRB_H
#define SCBRB_H
#include"Snap.h"
#include"utility.h"
#include"HeuMethods.h"

extern PUNGraph G;		//total graph
extern PUNGraph H;		//the current best community
extern int k_hat;		//the current biggest&feasible k
extern int low, high;	//the size bound

void BRB(TIntV C, TIntV R);
void Solve(int query, int l, int h);
int GetMinDegree(TIntV C);

#endif