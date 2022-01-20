#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <functional>
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> D2Array;
typedef IloArray<D2Array> D3Array;
typedef IloArray<D3Array> D4Array;
typedef IloArray<D4Array> D5Array;

FILE * covInfoFile = 0;
FILE * faultInfoFile = 0;
FILE * costInfoFile = 0;

int nTests = 0;
int nLines = 0;
int nFaults = 0;

int *lineFlag = 0;
int *faultFlag = 0;
int *faultID = 0;

int nTopTen = 0;
int **lineCover = 0;
int *topTenLines = 0;
int *topTenCover = 0;

struct test {
	int id;
	char code[10];
	int time;
	int nLinecover;
	int *lines;
	int nFaultcover;
	int *faults;
} *test_list;

void exportLPModel4FTCMinNTests(int budget, int nFaultsCovered) {

	int i, j, k, l, m;
	int temp, temp2, temp3;
	int * selectedTests;
	int nSelectedTests;
	int objF;

	IloEnv env; //Intialize the Cplex Environment 
	try {
		IloModel model(env);

		//---------------------Decision variable Definition Starts here---------------------------------------------
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		IloNumVarArray l(env, nLines, 0, 1, ILOINT);
		//---------------------Decision variable Definition Ends here-----------------------------------------------

		//---------------------Constraint Definition Starts here----------------------------------------------------	
		IloRangeArray cons(env);
		for (j = 0; j < nLines; j++) {
			if (lineFlag[j] > 0) {
				IloExpr v1(env);
				for (i = 0; i < nTests; i++) {
					temp = 0;
					for (k = 0; k < test_list[i].nLinecover; k++) {
						if (test_list[i].lines[k] == j + 1) temp = 1;
					}
					if (temp) v1 += x[i];
				}
				v1 += l[j];
				cons.add(v1 >= 1);
				v1.end();
			}
		}
		model.add(cons);

		for (j = 0; j < nFaults; j++) {
			IloExpr v2(env);
			for (i = 0; i < nTests; i++) {
				temp = 0;
				for (k = 0; k < test_list[i].nFaultcover; k++) {
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2 += x[i];
			}
			v2 += f[j];
			cons.add(v2 >= 1);
			v2.end();
		}
		model.add(cons);

		IloExpr v3(env);
		for (i = 0; i < nTests; i++) {
			v3 += x[i];
		}
		cons.add(v3 <= budget);
		v3.end();
		model.add(cons);

		IloExpr v4(env);
		for (i = 0; i < nLines; i++) {
			if (lineFlag[j] > 0) {
				v4 += l[i];
			}
		}
		cons.add(v4 == 0);
		v4.end();
		model.add(cons);

		IloExpr v5(env);
		for (i = 0; i < nFaults; i++) {
			v5 += f[i];
		}
		cons.add(v5 == nFaults - nFaultsCovered);
		v5.end();
		model.add(cons);
		//---------------------Constraint Definition Ends here-------------------------------------------------------

		//---------------------Objective Function--------------------------------------------------------------------
		IloExpr obj(env);
		for (j = 0; j < nTests; j++) {
			obj += x[j];
		}
		model.add(IloMinimize(env, obj));
		//---------------------Objective Function Def ends here----------------------------------------------------------

		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 21600);
		cplex.exportModel("lpexFTCminNTests.lp");
	} catch (...) {
		cerr << "Error" << endl;
	}
	env.end();
}

void exportLPModel4FTC(int budget) {

	int i, j, k, l, m;
	int temp, temp2, temp3;
	int * selectedTests;
	int nSelectedTests;
	int objF;

	IloEnv env; //Intialize the Cplex Environment 
	try {
		IloModel model(env);

		//---------------------Decision variable Definition Starts here---------------------------------------------
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		IloNumVarArray l(env, nLines, 0, 1, ILOINT);
		//---------------------Decision variable Definition Ends here-----------------------------------------------

		//---------------------Constraint Definition Starts here----------------------------------------------------	
		IloRangeArray cons(env);
		for (j = 0; j < nLines; j++) {
			if (lineFlag[j] > 0) {
				IloExpr v1(env);
				for (i = 0; i < nTests; i++) {
					temp = 0;
					for (k = 0; k < test_list[i].nLinecover; k++) {
						if (test_list[i].lines[k] == j + 1) temp = 1;
					}
					if (temp) v1 += x[i];
				}
				v1 += l[j];
				cons.add(v1 >= 1);
				v1.end();
			}
		}
		model.add(cons);

		for (j = 0; j < nFaults; j++) {
			IloExpr v2(env);
			for (i = 0; i < nTests; i++) {
				temp = 0;
				for (k = 0; k < test_list[i].nFaultcover; k++) {
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2 += x[i];
			}
			v2 += f[j];
			cons.add(v2 >= 1);
			v2.end();
		}
		model.add(cons);

		IloExpr v3(env);
		for (i = 0; i < nTests; i++) {
			v3 += x[i];
		}
		cons.add(v3 <= budget);
		v3.end();
		model.add(cons);
		//---------------------Constraint Definition Ends here-------------------------------------------------------

		//---------------------Objective Function--------------------------------------------------------------------
		IloExpr obj(env);
		for (j = 0; j < nFaults; j++) {
			obj += f[j];
		}
		for (j = 0; j < nLines; j++) {
			obj += l[j];
		}
		model.add(IloMinimize(env, obj));
		//---------------------Objective Function Def ends here----------------------------------------------------------
 
		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 21600);
		cplex.exportModel("lpexFTC.lp");
	} catch (...) {
		cerr << "Error" << endl;
	}
	env.end();
}

void exportLPModel4FVB() {

	int i, j, k, l, m;
	int temp, temp2, temp3;
	int * selectedTests;
	int nSelectedTests;

	IloEnv env; //Intialize the Cplex Environment 
	try {
		IloModel model(env);

		//---------------------Decision variable Definition Starts here---------------------------------------------
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		//---------------------Decision variable Definition Ends here-----------------------------------------------

		//---------------------Constraint Definition Starts here----------------------------------------------------	
		IloRangeArray cons(env);
		for (j = 0; j < nLines; j++) {
			if (lineFlag[j] > 0) {
				IloExpr v1(env);
				for (i = 0; i < nTests; i++) {
					temp = 0;
					for (k = 0; k < test_list[i].nLinecover; k++) {
						if (test_list[i].lines[k] == j + 1) temp = 1;
					}
					if (temp) v1 += x[i];
				}

				cons.add(v1 >= 1);
				v1.end();
			}
		}
		model.add(cons);

		for (j = 0; j < nFaults; j++) {
			IloExpr v2(env);
			for (i = 0; i < nTests; i++) {
				temp = 0;
				for (k = 0; k < test_list[i].nFaultcover; k++) {
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2 += x[i];
			}
			v2 += f[j];
			cons.add(v2 >= 1);
			v2.end();
		}
		model.add(cons);

		for (j = 0; j < nTopTen; j++) {
			if (lineFlag[topTenLines[j]] > 0) {
				IloExpr v3(env);
				for (i = 0; i < nTests; i++) {
					temp = 0;
					for (k = 0; k < test_list[i].nLinecover; k++) {
						if (test_list[i].lines[k] == topTenLines[j] + 1) temp = 1;
					}
					if (temp) v3 += x[i];
				}
				cons.add(v3 >= topTenCover[j]);
				v3.end();
			}
		}
		model.add(cons);
		//---------------------Constraint Definition Ends here-------------------------------------------------------

		//---------------------Objective Function--------------------------------------------------------------------
		IloExpr obj(env);
		for (j = 0; j < nFaults; j++) {
			obj += f[j];
		}
		for (j = 0; j < nTests; j++) {
			obj += x[j];
		}
		model.add(IloMinimize(env, obj));
		//---------------------Objective Function Def ends here----------------------------------------------------------
  	
		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 21600);
		cplex.exportModel("lpexFVB.lp");
	} catch (...) {
		cerr << "Error" << endl;
	}
	env.end();
}

void exportLPModel4FCB() {

	int i, j, k, l, m;
	int temp, temp2, temp3;
	int * selectedTests;
	int nSelectedTests;

	IloEnv env; //Intialize the Cplex Environment 
	try {
		IloModel model(env);

		//---------------------Decision variable Definition Starts here---------------------------------------------
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		//---------------------Decision variable Definition Ends here-----------------------------------------------

		//---------------------Constraint Definition Starts here----------------------------------------------------	
		IloRangeArray cons(env);
		for (j = 0; j < nLines; j++) {
			if (lineFlag[j] > 0) {
				IloExpr v1(env);
				for (i = 0; i < nTests; i++) {
					temp = 0;
					for (k = 0; k < test_list[i].nLinecover; k++) {
						if (test_list[i].lines[k] == j + 1) temp = 1;
					}
					if (temp) v1 += x[i];
				}

				cons.add(v1 >= 1);
				v1.end();
			}
		}
		model.add(cons);

		for (j = 0; j < nFaults; j++) {
			IloExpr v2(env);
			for (i = 0; i < nTests; i++) {
				temp = 0;
				for (k = 0; k < test_list[i].nFaultcover; k++) {
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2 += x[i];
			}
			v2 += f[j];
			cons.add(v2 >= 1);
			v2.end();
		}
		model.add(cons);
		//---------------------Constraint Definition Ends here-------------------------------------------------------

		//---------------------Objective Function--------------------------------------------------------------------
		IloExpr obj(env);
		for (j = 0; j < nFaults; j++) {
			obj += f[j];
		}
		for (j = 0; j < nTests; j++) {
			obj += x[j];
		}
		model.add(IloMinimize(env, obj));
		//---------------------Objective Function Def ends here----------------------------------------------------------

		IloCplex cplex(model);
		cplex.setParam(IloCplex::TiLim, 21600);
		cplex.exportModel("lpexFCB.lp");
	} catch (...) {
		cerr << "Error" << endl;
	}
	env.end();
}

char *trim(char *str) {
	size_t len = 0;
	char *frontp = str;
	char *endp = NULL;

	if (str == NULL) { return NULL; }
	if (str[0] == '\0') { return str; }

	len = strlen(str);
	endp = str + len;

	while (isspace((unsigned char)*frontp)) { ++frontp; }
	if (endp != frontp) {
		while (isspace((unsigned char) *(--endp)) && endp != frontp) {}
	}

	if (str + len - 1 != endp)
		*(endp + 1) = '\0';
	else if (frontp != str &&  endp == frontp)
		*str = '\0';

	endp = str;
	if (frontp != str) {
		while (*frontp) { *endp++ = *frontp++; }
		*endp = '\0';
	}
	return str;
}

int compare_ints(const void *a, const void *b) {
	if ((*(int **)a)[1] < (*(int **)b)[1]) return 1;
	else return -1;
}

void readdata() {
	int i, j, k;
	int tempInt, tempInt1, tempInt2;
	int lines;
	int ch = 0;
	int * templist;
	const char EOL = '\n';
	char str[15];
	int tempFault;
	int faultFl;
	int nLinesCovered;

	nTests = 0;
	while ((ch = fgetc(costInfoFile)) != EOF) 
		if (ch == '\n')
			nTests++;
	
	fclose(costInfoFile);
	test_list = (struct test *) malloc(nTests * sizeof(struct test));

	costInfoFile = fopen("rtime.info", "r");
	for (i = 0; i < nTests; i++) {
		test_list[i].id = i;
		fscanf(costInfoFile, "%s", &test_list[i].code);
		fscanf(costInfoFile, "%d", &test_list[i].time);
	}
	fclose(costInfoFile);

	tempInt = 0;
	while ((ch = fgetc(covInfoFile)) != EOF) 
		if (ch == '\n')
			tempInt++;
	
	fclose(covInfoFile);

	covInfoFile = fopen("cov.info", "r");
	for (i = 0; i < nTests; i++) {
		fscanf(covInfoFile, "%s", &str);
		tempInt = -1;
		for (j = 0; j < nTests; j++) 
			if (strcmp(test_list[j].code, str) == 0) 
				tempInt = j;
		
		tempInt2 = 0;
		while ((ch = fgetc(covInfoFile)) != EOL) 
			if (ch == ' ')
				tempInt2++;
		
		test_list[tempInt].nLinecover = tempInt2;
		test_list[tempInt].lines = (int *)malloc(test_list[tempInt].nLinecover * sizeof(int));
	}
	fclose(covInfoFile);

	covInfoFile = fopen("cov.info", "r");
	for (i = 0; i < nTests; i++) {
		fscanf(covInfoFile, "%s", &str);
		tempInt = -1;
		for (j = 0; j < nTests; j++)
			if (strcmp(test_list[j].code, str) == 0) 
				tempInt = j;

		for (j = 0; j < test_list[tempInt].nLinecover; j++) 
			fscanf(covInfoFile, "%d", &test_list[tempInt].lines[j]);
	}
	fclose(covInfoFile);
	lineFlag = (int *) calloc((size_t)(nLines), (size_t)sizeof(int));

	for (i = 0; i < nTests; i++)
		for (j = 0; j < test_list[i].nLinecover; j++)
			lineFlag[test_list[i].lines[j] - 1] = 1;

	tempInt = 0;
	while ((ch = fgetc(faultInfoFile)) != EOF)
		if (ch == '\n')
			tempInt++;
	
	fclose(faultInfoFile);
	templist = (int *)calloc((size_t)(nTests), (size_t)sizeof(int));

	faultInfoFile = fopen("fault.info", "r");
	for (i = 0; i < nTests; i++) {
		char buffer[20000];
		int numbers[10000];
		char *ptr;
		int cnt = 0;

		fgets(&buffer[0], 20000, faultInfoFile); // Get the string
		printf("%s\n", trim(buffer));
		ptr = strtok(buffer, " "); // Split in every space
		do {
			numbers[cnt] = strtol(ptr, NULL, 0); // Grab number
			if (strcmp(ptr, " ") != 0) 
				cnt++;
		} while ((ptr = strtok(NULL, " ")));
		templist[i] = cnt - 1;
	}
	fclose(faultInfoFile);
	faultID = (int *)calloc((size_t)(nFaults), (size_t)sizeof(int));

	faultInfoFile = fopen("fault.info", "r");
	tempFault = 0;
	for (i = 0; i < nTests; i++) {
		fscanf(faultInfoFile, "%s", &str);
		tempInt = -1;
		for (j = 0; j < nTests; j++)
			if (strcmp(test_list[j].code, str) == 0) 
				tempInt = j;

		test_list[tempInt].nFaultcover = templist[i];
		test_list[tempInt].faults = (int *)malloc(test_list[tempInt].nFaultcover * sizeof(int));

		for (j = 0; j < test_list[tempInt].nFaultcover; j++) {
			fscanf(faultInfoFile, "%d", &test_list[tempInt].faults[j]);
			faultFl = 0;
			for (k = 0; k < tempFault; k++)
				if (faultID[k] == test_list[tempInt].faults[j])
					faultFl = 1;
			
			if (faultFl < 1) {
				faultID[tempFault] = test_list[tempInt].faults[j];
				tempFault += 1;
			}
		}
	}

	fclose(covInfoFile);
	nTopTen = nLinesCovered*0.1;
	topTenLines = (int *)calloc((size_t)(nTopTen), (size_t)sizeof(int));
	topTenCover = (int *)calloc((size_t)(nTopTen), (size_t)sizeof(int));
	lineCover = (int *)calloc((size_t)(nLines), (size_t)sizeof(int));
	for (i = 0; i < nLines; i++) 
		lineCover[i] = (int *)calloc((size_t)(2), (size_t)sizeof(int));

	for (i = 0; i < nTests; i++) 
		for (j = 0; j < test_list[i].nLinecover; j++) 
			lineCover[test_list[i].lines[j] - 1][1] += 1;

	for (i = 0; i < nLines; i++)
		lineCover[i][0] = i;

	qsort(lineCover, nLines, sizeof(int*), compare_ints);

	for (i = 0; i < nTopTen; i++) {
		topTenLines[i] = lineCover[i][0];
		topTenCover[i] = 0.1*lineCover[i][1];
	}
}

void main(int argc, char *argv[]) {

	nLines = atoi(argv[1]); ;
	nFaults = atoi(argv[2]); ;

	covInfoFile = fopen("cov.info", "r");
	faultInfoFile = fopen("fault.info", "r");
	costInfoFile = fopen("rtime.info", "r");
	readdata();

	exportLPModel4FCB();
	exportLPModel4FVB();

	int budget = nTests*0.05;
	double tempDouble = (double)(nTests*0.05);
	if (tempDouble - budget > 0.5) budget += 1;
	exportLPModel4FTC(budget); // budget is set as 5% of all tests in this sample usage of the function

	exportLPModel4FTCMinNTests(budget, 30); // 30 is the # of faults to be covered as a constraint
}


