/*
 *  functions.h
 *  Distribution_1stEch
 *
 *  Created by Guven Ince on 3/31/13.
 *  Copyright 2013 UMASS. All rights reserved.
 *
 */

#ifndef DISASTER_ALLOCATION_FXNS_H_
#define DISASTER_ALLOCATION_FXNS_H_

#include "definitions.h"

vector<DNode> getDNodes();
vector<SNode> getSNodes();
vector<Cluster> getClusters();

void toLatexDemand(IloNumArray, vector<DNode>);
//void toLatexMIndex(IloNumArray2, vector<DNode>);
//void toLatexDuration(IloNumArray2, vector<SNode>, vector<DNode>);
void toLatexX(NumVar3dMatrix, vector<SNode>, vector<DNode>, IloCplex);
void toLatexP(IloCplex cplex, vector<DNode>);
//void toLatexYf(IntVar3dMatrix, vector<SNode>, vector<DNode>, IloCplex);
//void toLatexYb(IntVar3dMatrix, vector<SNode>, vector<DNode>, IloCplex);
void toLatexI(NumVarMatrix, vector<SNode>, IloCplex);
//void toLatexI1(IntVarMatrix, vector<SNode>, IloCplex);
//void toLatexI2(IntVarMatrix, vector<DNode>, IloCplex);
void toLatexU(NumVarMatrix, vector<DNode>, IloCplex);
void toLatexM(NumVarMatrix, vector<DNode>, IloCplex);
//void toLatexV(IntVarMatrix, vector<SNode>, IloCplex);
void toLatexS(NumVarMatrix, vector<SNode>, IloCplex);
//void toLatexStats(IloExpr, IloExpr, IloExpr, IloCplex);

#endif /* DISASTER_ALLOCATION_FXNS_H_ */
