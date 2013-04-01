// 6'dan farki 0 seviyesine geliste unutulan bir constraint eklendi
// Bu type1 ve type2 urunleri dagitabilen, Type1 urunleri replenishment sonrasi 0'da birakan final model.
// Version 2'den farki burada Agha'nin israrla koydurttugu bir constraint ve variable var
// Versiyon 3'dan farki bu model Type 1 urunleri de dagitabilir. Ancak bu urunler icin Tau'nun sifir olmasi gerekiyor
// Versiyon 4'un devami -- Inventory eklendi
#include <ilcplex/ilocplex.h>
#include <set>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

ILOSTLBEGIN

#include "functions.h"

using namespace std;

int main() {
#include "definitions.h"
	try {
		
		const char* datafile  = "data.txt";
		ifstream dfile(datafile);
		if ( !dfile ) {
			cerr << "ERROR: could not open file '" << datafile
			<< "' for reading" << endl;
			throw (-1);
		}
		
		dfile >> P >> S >> teta >> tau >> pComm >> goal >> pGoal;
		
		vector<DNode> districts = getDNodes(); // AVC, BAH, BAK, ..., SIL;
		vector<SNode> suppliers = getSNodes();
		vector<Cluster> clusters = getClusters();
		
		for (t=0; t < T; t++) {
			I[t] = NumVarMatrix(env, M);
			s[t] = NumVarMatrix(env, M);
			sT[t] = IloNumVarArray(env, K, 0, +IloInfinity);
			J[t] = IntVarMatrix(env, N);
			for(i=0; i < M; i++) {
				s[t][i] = IloNumVarArray(env, K, 0, +IloInfinity);
				I[t][i] = IloNumVarArray(env, K, 0, +IloInfinity);
			}
			for(j=0; j < N; j++)
				J[t][j] = IloIntVarArray(env, B, 0, 5*P[j]);
		}
		IloNumVarArray Sk(env, K, 0, +IloInfinity);
		
		for (t=1; t < T; t++) {
			p[t] = IntVar4dMatrix(env, N);
			m[t] = IntVarMatrix(env, N);
			for(j=0; j < N; j++) {
				p[t][j] = IntVar3dMatrix(env, teta[0]+tau[0]+1);
				m[t][j] = IloIntVarArray(env, K, 0, P[j]);
				for(a1 = 0; a1 <= teta[0]+tau[0]; a1++) {
					p[t][j][a1] = IntVarMatrix(env, teta[1]+tau[1]+1);
					for(a2 = 1; a2 <= teta[1]+tau[1]; a2++){
						p[t][j][a1][a2] = IloIntVarArray(env, teta[2]+tau[2]+1);
						for(a3 = 1; a3 <= teta[2]+tau[2]; a3++){
							p[t][j][a1][a2][a3] = IloIntVar(env, 0, P[j]);
						}
					}
				}
			}
		}
		
		for (t=1; t < T; t++) {
			x[t] = NumVar3dMatrix(env, M);
			X[t] = IntVar5dMatrix(env, N);
			E[t] = IntVarMatrix(env, N);
			zeta[t] = IntVar3dMatrix(env, N);
			for(j=0; j < N; j++) {
				X[t][j] = IntVar4dMatrix(env, teta[0]+tau[0]);
				E[t][j] =  IloIntVarArray(env, B, 0, P[j]);
				zeta[t][j] = IntVarMatrix(env, K);
				for(a1 = 0; a1 < teta[0]+tau[0]; a1++) { //0'dan mi baslatmak lazim acaba
					X[t][j][a1] = IntVar3dMatrix(env, teta[1]+tau[1]);
					for(a2 = 1; a2 < teta[1]+tau[1]; a2++){
						X[t][j][a1][a2] = IntVarMatrix(env, teta[2]+tau[2]);
						for(a3 = 1; a3 < teta[2]+tau[2]; a3++)
							X[t][j][a1][a2][a3] = IloIntVarArray(env, B, 0, P[j]);
					}
				}
				for(k=0; k < K; k++) {
					zeta[t][j][k] = IloIntVarArray(env, B, 0, P[j]);
				}
			}
			for(i=0; i < M; i++) {
				x[t][i] = NumVarMatrix(env, N);
				for(j=0; j < N; j++) {
					x[t][i][j] = IloNumVarArray(env, K, 0, P[j]);
				}
			}
		}
		
		//OBJECTIVE0 :: Minimize total mortality
		IloExpr tMort(env);
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (k = 0; k < K; k++)
					tMort += m[t][j][k];
		
		//OBJECTIVE1 :: Minimize total suffering population
		IloExpr tSuff1(env);
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = tau[0] + 1; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 < teta[2]+tau[2]; a3++)
							tSuff1 += p[t][j][a1][a2][a3] - X[t][j][a1][a2][a3][1] - X[t][j][a1][a2][a3][4] - X[t][j][a1][a2][a3][5] - X[t][j][a1][a2][a3][7];
		IloExpr tSuff2(env);
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
					for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 < teta[2]+tau[2]; a3++)
							tSuff2 += p[t][j][a1][a2][a3] - X[t][j][a1][a2][a3][2] - X[t][j][a1][a2][a3][4] - X[t][j][a1][a2][a3][6] - X[t][j][a1][a2][a3][7];
		IloExpr tSuff3(env);
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++)
							tSuff3 += p[t][j][a1][a2][a3] - X[t][j][a1][a2][a3][3] - X[t][j][a1][a2][a3][5] - X[t][j][a1][a2][a3][6] - X[t][j][a1][a2][a3][7];;
		IloExpr tSuff(env);
		tSuff = pComm[0]*tSuff1 + pComm[1]*tSuff2 + pComm[2]*tSuff3;
		
		//OBJECTIVE2 :: Minimize total supply
		IloExpr tSk(env);
		for (k = 0; k < K; k++)
			tSk += Sk[k];
		
		//GOAL OBJECTIVE
		IloExpr tGoal(env);
		tGoal = pGoal[0]*d1 + pGoal[1]*d2;
		
		vector<IloObjective> objs;
		vector<IloExpr> Z;
		Z.push_back(tMort);
		Z.push_back(tSuff);
		Z.push_back(tSk);
		Z.push_back(tGoal);
		
		IloObjective minMort(env, Z[0], IloObjective::Minimize);
		IloObjective minSuff(env, Z[1], IloObjective::Minimize);
		IloObjective minSk(env, Z[2], IloObjective::Minimize);
		IloObjective minGoal(env, Z[3], IloObjective::Minimize);
		objs.push_back(minMort);
		objs.push_back(minSuff);
		objs.push_back(minSk);
		objs.push_back(minGoal);
		
		cout << "RESTRICTIONS" << endl;
		for(t = 1; t < T; t++)
			for (j = 0; j < N; j++){
				for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 <= tau[1]; a2++)
						for (a3 = 1; a3 < teta[2]+tau[2]; a3++)
							mod.add(X[t][j][a1][a2][a3][2] + X[t][j][a1][a2][a3][4] + X[t][j][a1][a2][a3][6] + X[t][j][a1][a2][a3][7] == 0);
				for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 <= tau[2]; a3++)
							mod.add(X[t][j][a1][a2][a3][3] + X[t][j][a1][a2][a3][5] + X[t][j][a1][a2][a3][6] + X[t][j][a1][a2][a3][7] == 0);
				for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
					for (a3 = 1; a3 < teta[2]+tau[2]; a3++)
						mod.add(X[t][j][0][a2][a3][1] + X[t][j][0][a2][a3][4] + X[t][j][0][a2][a3][5] + X[t][j][0][a2][a3][7] == 0);
			}
		/*
		 for(k=0; k<K;k++) {
		 mod.add(sT[0][k] == S[k]/2);
		 }
		 int ts1 = 12;
		 int ts2 = 10;
		 int ts3	= 6;
		 for(t=1; t<=ts1;t++) {
		 mod.add(sT[t][0] == t*(S[0]/(ts1*(ts1+1))));
		 }
		 for(t=1; t<=ts2;t++) {
		 mod.add(sT[t][1] == t*(S[1]/(ts2*(ts2+1))));
		 }
		 for(t=1; t<=ts3;t++) {
		 mod.add(sT[t][2] == t*(S[2]/(ts3*(ts3+1))));
		 }
		 */
		//Adding objectives to env
		mod.add(objs[1]);
		
		cout << "GOAL CONSTRAINTS " << endl;
		//				mod.add(tMort == 0);
		//				mod.add(tSuff == 0);
		
		cout << "CONSTRAINT 1" << endl;
		for (k = 0; k < K; k++) {
			IloExpr tSply(env);
			for (t = 0; t < T; t++)
				tSply += sT[t][k];
			mod.add(tSply <= S[k]);
			tSply.end();
		}
		
		cout << "CONSTRAINT 2" << endl;
		for (t = 0; t < T; t++) {
			for (k = 0; k < K; k++) {
				IloExpr sOi(env);
				for (i = 0; i < M; i++)
					sOi += s[t][i][k];
				mod.add(sOi - sT[t][k] == 0);
				sOi.end();
			}
		}
		
		cout << "CONSTRAINT 3" << endl;
		//initialization of initial inventory to 0
		for (k = 0; k < K; k++)
			for (i = 0; i < M; i++)
				mod.add(I[0][i][k] - s[0][i][k]== 0);
		
		cout << "CONSTRAINT 4" << endl;
		// Flow balance around supply points
		for (t = 1; t < T; t++) {
			for (k = 0; k < K; k++) {
				for (i = 0; i < M; i++) {
					IloExpr xOj(env);
					for (j = 0; j < N; j++) {
						xOj += x[t][i][j][k];
					}
					mod.add(I[t-1][i][k] + s[t][i][k] - xOj - I[t][i][k] == 0);
					xOj.end();
				}
			}
		}
		
		cout << "CONSTRAINT 5" << endl;
		//agha'nin gereksiz konstrainti. 6 ile birlesebilir, Versiyon 8'de oldugu gibi.
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (k = 0; k < K; k++){
					IloExpr XOi(env);
					for (i = 0; i < M; i++)
						XOi += x[t][i][j][k];
					if (k == 0)
						mod.add(XOi == zeta[t][j][k][1] + zeta[t][j][k][4] + zeta[t][j][k][5] + zeta[t][j][k][7]);
					else if (k == 1)
						mod.add(XOi == zeta[t][j][k][2] + zeta[t][j][k][4] + zeta[t][j][k][6] + zeta[t][j][k][7]);
					else
						mod.add(XOi == zeta[t][j][k][3] + zeta[t][j][k][5] + zeta[t][j][k][6] + zeta[t][j][k][7]);
					XOi.end();
				}
		
		cout << "CONSTRAINT 6" << endl;
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (k = 0; k < K; k++){
					if (k == 0){
						mod.add(zeta[t][j][k][1] == E[t][j][1]);
						mod.add(zeta[t][j][k][4] == E[t][j][4]);
						mod.add(zeta[t][j][k][5] == E[t][j][5]);
						mod.add(zeta[t][j][k][7] == E[t][j][7]);
					}
					else if (k == 1){
						mod.add(zeta[t][j][k][2] == E[t][j][2]);
						mod.add(zeta[t][j][k][4] == E[t][j][4]);
						mod.add(zeta[t][j][k][6] == E[t][j][6]);
						mod.add(zeta[t][j][k][7] == E[t][j][7]);
					}
					else if (k == 2){
						mod.add(zeta[t][j][k][3] == E[t][j][3]);
						mod.add(zeta[t][j][k][5] == E[t][j][5]);
						mod.add(zeta[t][j][k][6] == E[t][j][6]);
						mod.add(zeta[t][j][k][7] == E[t][j][7]);
					}
				}
		
		cout << "CONSTRAINT 7-ek" << endl;
		//initialization of initial inventory to 0
		for (b = 1; b < B; b++)
			for (j = 0; j < N; j++)
				mod.add(J[0][j][b] == 0);
		
		cout << "CONSTRAINT 7" << endl;
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for(b=1; b < B; b++){
					IloExpr XOa123(env);
					for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
						for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
							for (a3 = 1; a3 < teta[2]+tau[2]; a3++)
								XOa123 += X[t][j][a1][a2][a3][b];
					mod.add(J[t-1][j][b] + E[t][j][b] - XOa123 - J[t][j][b] == 0);
					XOa123.end();
				}
		
		cout << "CONSTRAINT 8" << endl;
		//initialization to population
		for (j = 0; j < N; j++)
			for (a1 = 0; a1 <= teta[0]+tau[0]; a1++)
				for (a2 = 1; a2 <= teta[1]+tau[1]; a2++)
					for (a3 = 1; a3 <= teta[2]+tau[2]; a3++){
						if ((a1 == tau[0]+1) && (a2 == tau[1]+1) && (a3 == tau[2]+1))
							mod.add(p[1][j][a1][a2][a3] == P[j]);
						else
							mod.add(p[1][j][a1][a2][a3] == 0);
					}
		
		cout << "CONSTRAINT 9" << endl;
		//if bundle includes only 1 group goes to 0,a2a3 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a2 = 2; a2 <= teta[1]+tau[1]; a2++)
					for (a3 = 2; a3 <= teta[2]+tau[2]; a3++){
						IloExpr xOa1(env);
						for (a1 = tau[0]+1; a1 < teta[0]+tau[0]; a1++){
							xOa1 += X[t-1][j][a1][a2-1][a3-1][1];
						}
						IloExpr xOa10(env);
						for (b = 1; b < B; b++){
							xOa10 += X[t-1][j][0][a2-1][a3-1][b];
						}
						mod.add(p[t][j][0][a2][a3] == xOa1 + p[t-1][j][0][a2-1][a3-1] - xOa10);
						xOa1.end();
						xOa10.end();
					}
		
		cout << "CONSTRAINT 10" << endl;
		//if bundle includes only 2 group goes to a1,1,a3 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 2; a1 <= teta[0]+tau[0]; a1++)
					for (a3 = 2; a3 <= teta[2]+tau[2]; a3++){
						IloExpr xOa2(env);
						for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++){
							xOa2 += X[t-1][j][a1-1][a2][a3-1][2];
						}
						mod.add(p[t][j][a1][1][a3] == xOa2);
						xOa2.end();
					}
		
		cout << "CONSTRAINT 11" << endl;
		//if bundle includes only 3 group goes to a1a2,1 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 2; a1 <= teta[0]+tau[0]; a1++)
					for (a2 = 2; a2 <= teta[1]+tau[1]; a2++){
						IloExpr xOa3(env);
						for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++){
							xOa3 += X[t-1][j][a1-1][a2-1][a3][3];
						}
						mod.add(p[t][j][a1][a2][1] == xOa3);
						xOa3.end();
					}
		
		cout << "CONSTRAINT 12" << endl;
		//if bundle includes 1 & 2 group goes to 01,a3 level for a1 >1; if bundle includes only 2 for a1=1
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a3 = 2; a3 <= teta[2]+tau[2]; a3++) {
					IloExpr xOa12(env);
					for (a1 = tau[0] + 1; a1 < teta[0]+tau[0]; a1++)
						for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++){
							xOa12 += X[t-1][j][a1][a2][a3-1][4];
						}
					IloExpr xOa2(env);
					for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++){
						xOa2 += X[t-1][j][0][a2][a3-1][2];
					}
					mod.add(p[t][j][0][1][a3] == xOa12 + xOa2);
					xOa12.end();
					xOa2.end();
				}
		
		cout << "CONSTRAINT 13" << endl;
		//if bundle includes 1 & 3 group goes to 0,a2,1 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a2 = 2; a2 <= teta[1]+tau[1]; a2++) {
					IloExpr xOa13(env);
					for (a1 = tau[0] + 1; a1 < teta[0]+tau[0]; a1++)
						for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++){
							xOa13 += X[t-1][j][a1][a2-1][a3][5];
						}
					IloExpr xOa3(env);
					for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++){
						xOa3 += X[t-1][j][0][a2-1][a3][3];
					}
					mod.add(p[t][j][0][a2][1] == xOa13 + xOa3);
					xOa13.end();
					xOa3.end();
				}
		
		cout << "CONSTRAINT 14" << endl;
		//if bundle includes 2 & 3 group goes to a1,11 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 2; a1 <= teta[0]+tau[0]; a1++) {
					IloExpr xOa23(env);
					for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++){
							xOa23 += X[t-1][j][a1-1][a2][a3][6];
						}
					mod.add(p[t][j][a1][1][1] == xOa23);
					xOa23.end();
				}
		
		cout << "CONSTRAINT 15" << endl;
		//if bundle includes all commodities group goes to 011 level.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++){
				IloExpr xOa123(env);
				for (a1 = tau[0] + 1; a1 < teta[0]+tau[0]; a1++)
					for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++)
							xOa123 += X[t-1][j][a1][a2][a3][7];
				IloExpr xOa23(env);
				for (a2 = tau[1] + 1; a2 < teta[1]+tau[1]; a2++)
					for (a3 = tau[2] + 1; a3 < teta[2]+tau[2]; a3++)
						xOa23 += X[t-1][j][0][a2][a3][6];
				mod.add(p[t][j][0][1][1] == xOa123 + xOa23);
				xOa123.end();
				xOa23.end();
			}
		
		cout << "CONSTRAINT 16" << endl;
		//without shipment, group goes to upper level in all comm.
		for (t = 2; t < T; t++)
			for (j = 0; j < N; j++){
				for (a1 = 2; a1 <= tau[0]+teta[0]; a1++)
					for (a2 = 2; a2 <= tau[1]+teta[1]; a2++)
						for (a3 = 2; a3 <= tau[2]+teta[2]; a3++){
							IloExpr XOb(env);
							for (b = 1; b < B; b++)
								XOb += X[t-1][j][a1-1][a2-1][a3-1][b];
							mod.add(p[t][j][a1][a2][a3] == p[t-1][j][a1-1][a2-1][a3-1] - XOb);
						}
				/*
				 for (a2 = 2; a2 <= tau[1]+teta[1]; a2++)
				 mod.add(p[t][j][0][a2][tau[2]+1] == p[t-1][j][0][a2-1][tau[2]]);
				 for (a3 = 2; a3 <= tau[2]+teta[2]; a3++)
				 mod.add(p[t][j][0][tau[1]+1][a3] == p[t-1][j][0][tau[1]][a3-1]);
				 mod.add(p[t][j][0][tau[1]+1][tau[2]+1] == p[t-1][j][0][tau[1]][tau[2]]);
				 for (a1 = 2; a1 <= teta[0]+tau[0]; a1++)
				 for (a2 = 2; a2 <= tau[1]+teta[1]; a2++)
				 mod.add(p[t][j][a1][a2][tau[2]+1] == p[t-1][j][a1-1][a2-1][tau[2]]);
				 for (a1 = 2; a1 <= teta[0]+tau[0]; a1++)
				 for (a3 = 2; a3 <= tau[2]+teta[2]; a3++)
				 mod.add(p[t][j][a1][tau[1]+1][a3] == p[t-1][j][a1-1][tau[1]][a3-1]);
				 for (a1 = 2; a1 <= teta[0]+tau[0]; a1++)
				 mod.add(p[t][j][a1][tau[1]+1][tau[2]+1] == p[t-1][j][a1-1][tau[1]][tau[2]]);
				 for (a2 = 2; a2 <= teta[1]+tau[1]; a2++)
				 mod.add(p[t][j][0][a2][tau[2]+1] == p[t-1][j][0][a2-1][tau[2]]);
				 for (a3 = 2; a3 <= teta[2]+tau[2]; a3++)
				 mod.add(p[t][j][0][tau[1]+1][a3] == p[t-1][j][0][tau[1]][a3-1]);
				 */
			}
		
		cout << "CONSTRAINT 17" << endl;
		//capturing mortality
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++) {
				{
					IloExpr pOa3(env);
					for (a1 = 0; a1 <= teta[0]+tau[0]; a1++)
						for (a2 = 1; a2 <= teta[1]+tau[1]; a2++)
							pOa3 += p[t][j][a1][a2][teta[2]+tau[2]];
					mod.add(m[t][j][2] == pOa3);
					pOa3.end();
				}
				{
					IloExpr pOa2(env);
					for (a1 = 0; a1 <= teta[0]+tau[0]; a1++)
						for (a3 = 1; a3 <= teta[2]+tau[2]; a3++)
							pOa2 += p[t][j][a1][teta[1]+tau[1]][a3];
					mod.add(m[t][j][1] == pOa2);
					pOa2.end();
				}
				{
					IloExpr pOa1(env);
					for (a2 = 1; a2 <= teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 <= teta[2]+tau[2]; a3++)
							pOa1 += p[t][j][teta[0]+tau[0]][a2][a3];
					mod.add(m[t][j][0] == pOa1);
					pOa1.end();
				}
			}
		
		cout << "CONSTRAINT 18" << endl;
		//if zero population, no shipment
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 0; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 < teta[2]+tau[2]; a3++){
							//							for (b = 1; b < B; b++)
							//								mod.add(X[t][j][a1][a2][a3][b] <= p[t][j][a1][a2][a3]);
							mod.add(X[t][j][a1][a2][a3][1] + X[t][j][a1][a2][a3][4] + X[t][j][a1][a2][a3][5] + X[t][j][a1][a2][a3][7] <= p[t][j][a1][a2][a3]);
							mod.add(X[t][j][a1][a2][a3][2] + X[t][j][a1][a2][a3][4] + X[t][j][a1][a2][a3][6] + X[t][j][a1][a2][a3][7] <= p[t][j][a1][a2][a3]);
							mod.add(X[t][j][a1][a2][a3][3] + X[t][j][a1][a2][a3][5] + X[t][j][a1][a2][a3][6] + X[t][j][a1][a2][a3][7] <= p[t][j][a1][a2][a3]);
						}
		
		cout << "CONSTRAINT 19" << endl; // bu constrainte dikkat!!!!!!!!!
																		 // a1 = 2 hicbir constraintte deger almiyor. bu yuzden sifirliyoruz.
		for(t = 2; t < T; t++)
			for (j = 0; j < N; j++)
				for (a2 = 1; a2 <= teta[1]+tau[1]; a2++)
					for (a3 = 1; a3 <= teta[2]+tau[2]; a3++){
						mod.add(p[t][j][1][a2][a3] == 0);
					}
		
		cout << "CONSTRAINT SON" << endl;
		cout << "DENEME OBJECTIVE" << endl;
		
		IloCplex cplex(env);
		cplex.extract(mod);
		cplex.exportModel("model.lp");
		
		cplex.solve();
		
		cplex.out() << "solution status = " << cplex.getStatus() << endl;
		cplex.out() << "objective value = " << cplex.getObjValue() << endl;
		
		cout << "TEST SON" << endl;
		
		for(j=0;j<N;j++)
			toLatexP(cplex, districts);
		vector<vector<int> > sufferers(3, vector<int>(13));
		int suff1 = 0; int t1 = 0;
		int suff2 = 0; int t2 = 0;
		int suff3 = 0; int t3 = 0;
		for (t = 1; t < T; t++)
			for (j = 0; j < N; j++)
				for (a1 = 1; a1 < teta[0]+tau[0]; a1++)
					for (a2 = 1; a2 < teta[1]+tau[1]; a2++)
						for (a3 = 1; a3 < teta[2]+tau[2]; a3++){
							if (cplex.getValue(p[t][j][a1][a2][a3]) > 0){
								if (a1 > tau[0]){
									suff1 += cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][1]) - cplex.getValue(X[t][j][a1][a2][a3][4]) - cplex.getValue(X[t][j][a1][a2][a3][5]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									sufferers[0][t1] = cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][1]) - cplex.getValue(X[t][j][a1][a2][a3][4]) - cplex.getValue(X[t][j][a1][a2][a3][5]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									t1 += 1;
								}
								if (a2 > tau[1]){
									suff2 += cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][2]) - cplex.getValue(X[t][j][a1][a2][a3][4]) - cplex.getValue(X[t][j][a1][a2][a3][6]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									sufferers[1][t2] = cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][2]) - cplex.getValue(X[t][j][a1][a2][a3][4]) - cplex.getValue(X[t][j][a1][a2][a3][6]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									t2 += 1;
								}
								if (a3 > tau[2]){
									suff3 += cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][3]) - cplex.getValue(X[t][j][a1][a2][a3][5]) - cplex.getValue(X[t][j][a1][a2][a3][6]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									sufferers[2][t3] = cplex.getValue(p[t][j][a1][a2][a3]) - cplex.getValue(X[t][j][a1][a2][a3][3]) - cplex.getValue(X[t][j][a1][a2][a3][5]) - cplex.getValue(X[t][j][a1][a2][a3][6]) - cplex.getValue(X[t][j][a1][a2][a3][7]);
									t3 += 1;
								}
							}
						}
		cout << "suff1 = " << suff1 << ", t1 = " << t1 << endl;
		cout << "suff2 = " << suff2 << ", t2 = " << t2 << endl;
		cout << "suff3 = " << suff3 << ", t3 = " << t3 << endl;
		cout << "---- sufferers -----" << endl;
		for (int l = 0; l < T; l++)
			cout << l << "\t"<< sufferers[0][l] << "\t" << sufferers[1][l] << "\t" << sufferers[2][l] << endl;
		
		cout << "sT[t][k]" << endl;
		for(t = 0; t < T; t++){
			cout << t << "\t";
			for(k=0; k<K; k++)
				cout <<  cplex.getValue(sT[t][k]) << "\t";
			cout << endl;
		}
		
		cout << "s[t][i][k]" << endl;
		for(t = 0; t < T; t++){
			cout << "t = " << t << endl;
			for(i = 0; i < M; i++){
				cout << "i = " << i << "\t";
				for(k=0; k<K; k++)
					cout <<  cplex.getValue(s[t][i][k]) << "\t";
				cout << endl;
			}
		}
		
		cout << "I[t][i][k]" << endl;
		for(t = 0; t < T; t++){
			cout << "t = " << t << endl;
			for(i = 0; i < M; i++){
				cout << "i = " << i << "\t";
				for(k=0; k<K; k++)
					cout <<  cplex.getValue(I[t][i][k]) << "\t";
				cout << endl;
			}
		}
		
		cout << "x[t][i][j][k]" << endl;
		for(t = 1; t < T; t++){
			cout << "t = " << t << endl;
			for(i = 0; i < M; i++){
				cout << "i = " << i << "\t";
				for(k=0; k<K; k++)
					for(j = 0; j < N; j++)
						cout <<  cplex.getValue(x[t][i][j][k]) << "\t";
				cout << endl;
			}
		}
		
		//		cout << "Sk[k]" << endl;
		//				for(k=0; k<K; k++)
		//					cout << fixed << setprecision(0) << cplex.getValue(Sk[k]) << "\t";
		//				cout << "total = " << cplex.getValue(tSk) << endl;
		
		cplex.out() << "objective value = " << cplex.getObjValue() << endl;
		cout << "tMort = " << cplex.getValue(tMort) << endl;
		cout << "tSuff = " << cplex.getValue(tSuff) << endl;
		
//		env.end();
//		mod.end();
	}
	catch (IloException& e) {
		cerr << "ERROR: " << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	return 0;
}
