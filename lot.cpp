#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ilcplex/ilocplex.h>
#include <iostream>

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <queue>
#include <algorithm>
#include <string>
#include <cassert>
#include <cmath>
#include <iomanip>

#include <functional>
#include <vector>
using std::vector;

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN


typedef IloArray<IloNumVarArray> D2Array;    //Define 2 Dimensional Array 
typedef IloArray<D2Array> D3Array;           //Define 3 Dimensional Array 
typedef IloArray<D3Array> D4Array;           //Define 4 Dimensional Array 
typedef IloArray<D4Array> D5Array;           //Define 4 Dimensional Array 

FILE * infile;
FILE * infile2;
FILE * infile3;
FILE * infile4;
FILE * outfile;
FILE * outfile2;

int nTests; 
int nLines; 
int nFaults; 

int *lineFlag; 
int *faultFlag; 
int *faultID; 

int nTopTen; 
int **lineCover;
int *topTenLines; 
int *topTenCover; 

int LC, FC; 

struct test
{
	int id;
	char code[10];
	int time;
	int nLinecover;
	int *lines;
	int nFaultcover;
	int *faults;

} *test_list;





double Solver4(int nTestCap){

	int i,j,k,l,m;
	int temp, temp2, temp3; 
	double tempDouble; 
	int * selectedTests; 
	int nSelectedTests; 
	int objF; 

		
	IloEnv env; //Intialize the Cplex Environment 
	try{
		IloModel model(env);
		
	//---------------------Decision variable Definition Starts here---------------------------------------------

		//IloNumVarArray Inv(env, planHorizon, 0, IloInfinity, ILOINT);
		//IloNumVarArray Prod(env, planHorizon, 0, IloInfinity, ILOINT);
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		IloNumVarArray l(env, nLines, 0, 1, ILOINT);

		//D2Array Inv(env,nItems); 
		//for(j=0;j<nItems;j++){
		//	Inv[j] = IloNumVarArray(env,planHorizon, 0, IloInfinity, ILOINT); 	
		//}
	
	//---------------------Decision variable Definition Ends here-----------------------------------------------

	//---------------------Constraint Definition Starts here----------------------------------------------------	
		
		IloRangeArray cons(env);
		for(j=0;j<nLines;j++){
			if (lineFlag[j]>0){
				IloExpr v1(env); 
				for(i=0;i<nTests;i++){
					temp = 0;
					for(k=0;k<test_list[i].nLinecover;k++){
						if (test_list[i].lines[k] == j+1) temp = 1;
					}
					if (temp) v1+=x[i];
				}
				v1+=l[j];	
				cons.add(v1>=1);
				v1.end();
			}
		}
		model.add(cons); 

		for(j=0;j<nFaults;j++){
			IloExpr v2(env); 
			for(i=0;i<nTests;i++){
				temp = 0;
				for(k=0;k<test_list[i].nFaultcover;k++){
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2+=x[i];
			}
			v2+=f[j];
			cons.add(v2>=1);
			v2.end();
		}
		model.add(cons); 

		IloExpr v3(env); 
		for(i=0;i<nTests;i++){
			v3+=x[i];
		}
		cons.add(v3<=nTestCap);
		v3.end();
		model.add(cons); 



		IloExpr v4(env); 
		for(i=0;i<nLines;i++){
			if (lineFlag[j]>0){
				v4+=l[i];
			}
		}
		cons.add(v4==0);
		v4.end();
		model.add(cons); 



		IloExpr v5(env); 
		for(i=0;i<nFaults;i++){
			v5+=f[i];
		}
		cons.add(v5==nFaults-FC);
		v5.end();
		model.add(cons); 
		
	//---------------------Constraint Definition Ends here-------------------------------------------------------

	//---------------------Objective Function--------------------------------------------------------------------

		IloExpr obj(env);

		//for(j=0;j<nFaults;j++){
		//	obj+=f[j];
		//}

		//for(j=0;j<nLines;j++){
		//	obj+=l[j];
		//}
		
		for(j=0;j<nTests;j++){
			obj+=x[j];
		}

		model.add(IloMinimize(env, obj));			
			
	//---------------------Objective Function Def ends here----------------------------------------------------------
		
	//---------------------Solving the model-------------------------------------------------------------------------  
			
		IloCplex cplex(model);
		//if (easyFlag == 1) cplex.setParam(IloCplex::EpGap, 0.8); 
		//else 
		//cplex.setParam(IloCplex::EpGap, 0.25); 
		//cplex.setParam(IloCplex::EpGap, 0.005); 
		//if (easyFlag == 1) cplex.setParam(IloCplex::TiLim, 10);
		//else cplex.setParam(IloCplex::TiLim, 120);
		//else 
		cplex.setParam(IloCplex::TiLim, 21600);
		//cplex.setOut(env.getNullStream());
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);    //Select barrier algorithm 
        //cplex.setParam(IloCplex::BarCrossAlg, IloCplex::NoAlg);  //Without crossover

		//cplex.exportModel("lpex.lp");
		cplex.solve();
		
	//----------------------Extracting results for output-------------------------------------------------------------	
		
	tempDouble = cplex.getObjValue();
	
	fprintf(outfile, "Obj. Fun. %lf\n",cplex.getObjValue());  
	fprintf(outfile, "Gap. %lf\n", cplex.getMIPRelativeGap());

	////fprintf(outfile, "%lf ",cplex.getObjValue());  

	fprintf(outfile, "Selected Tests\n"); 
	
	nSelectedTests = 0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		//if (temp>0) fprintf (outfile,"%d %s\n", j+1, test_list[j].code);  
		nSelectedTests += temp; 
	}
	
	fprintf(outfile, "Number of Selected Tests %d\n", nSelectedTests); 
	fprintf(outfile2, "%d ", nSelectedTests); 
	
	memset(selectedTests, 0, (unsigned)(nSelectedTests)* (unsigned)sizeof(int));
	temp2=0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		if (temp>0) {
			selectedTests[temp2] = j; 
			temp2 = temp2 + 1;  
		}
	}


	fprintf(outfile, "Checking Line Coverage\n"); 

	temp2=0; 
	for(j=0;j<nLines;j++){
		if (lineFlag[j]>0){
			temp = 0; 
			for(i=0;i<nSelectedTests;i++){
				//fprintf (outfile,"%d %s\n", selectedTests[i], test_list[selectedTests[i]].code);  
				for(k=0;k<test_list[selectedTests[i]].nLinecover;k++){
					if (test_list[selectedTests[i]].lines[k] ==j+1) temp =1; 
				}
			}
			if (temp <1) fprintf(outfile, "Line Not Covered %d\n", j+1); 
			temp2 = temp2 + temp; 
		}
	}

	fprintf(outfile, "Number of Lines Covered %d\n", temp2); 
	fprintf(outfile2, "%d ", temp2); 
	LC = temp2; 
	objF = temp2; 

	fprintf(outfile, "Faults\n"); 
	
	temp2 = 0; 
	for(j=0;j<nFaults;j++){
		temp = cplex.getValue(f[j]);
		if (cplex.getValue(f[j])-temp > 0.8) temp = temp+1;
		if (temp>0) fprintf (outfile,"%d\n", faultID[j]);  
		temp2 += temp; 
	}
	
	fprintf(outfile, "Uncovered Faults %d\n", temp2); 
	fprintf(outfile2, "%d ", nFaults-temp2); 
	FC = nFaults-temp2; 
	objF += nFaults-temp2; 
	fprintf(outfile2, "%d ", objF); 

	for(j=0;j<nFaults;j++){
		temp = 0;
		for(i=0;i<nSelectedTests;i++){
			for(k=0;k<test_list[selectedTests[i]].nFaultcover;k++){
				if (test_list[selectedTests[i]].faults[k] == faultID[j]) temp = 1;
				//fprintf(outfile, "%d %d\n", test_list[selectedTests[i]].faults[k], temp); 
			}
		}
		if (temp <1) fprintf(outfile, "Fault Not Covered %d\n", faultID[j]); 
	}
		


	}

   catch (IloException& e) {
      cerr << "Error: " << e.getMessage() << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

	env.end();

	return tempDouble;
	


}



double Solver3(int nTestCap){

	int i,j,k,l,m;
	int temp, temp2, temp3; 
	double tempDouble; 
	int * selectedTests; 
	int nSelectedTests; 
	int objF; 

		
	IloEnv env; //Intialize the Cplex Environment 
	try{
		IloModel model(env);
		
	//---------------------Decision variable Definition Starts here---------------------------------------------

		//IloNumVarArray Inv(env, planHorizon, 0, IloInfinity, ILOINT);
		//IloNumVarArray Prod(env, planHorizon, 0, IloInfinity, ILOINT);
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);
		IloNumVarArray l(env, nLines, 0, 1, ILOINT);

		//D2Array Inv(env,nItems); 
		//for(j=0;j<nItems;j++){
		//	Inv[j] = IloNumVarArray(env,planHorizon, 0, IloInfinity, ILOINT); 	
		//}
	
	//---------------------Decision variable Definition Ends here-----------------------------------------------

	//---------------------Constraint Definition Starts here----------------------------------------------------	
		
		IloRangeArray cons(env);
		for(j=0;j<nLines;j++){
			if (lineFlag[j]>0){
				IloExpr v1(env); 
				for(i=0;i<nTests;i++){
					temp = 0;
					for(k=0;k<test_list[i].nLinecover;k++){
						if (test_list[i].lines[k] == j+1) temp = 1;
					}
					if (temp) v1+=x[i];
				}
				v1+=l[j];	
				cons.add(v1>=1);
				v1.end();
			}
		}
		model.add(cons); 

		for(j=0;j<nFaults;j++){
			IloExpr v2(env); 
			for(i=0;i<nTests;i++){
				temp = 0;
				for(k=0;k<test_list[i].nFaultcover;k++){
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2+=x[i];
			}
			v2+=f[j];
			cons.add(v2>=1);
			v2.end();
		}
		model.add(cons); 

		IloExpr v3(env); 
		for(i=0;i<nTests;i++){
			v3+=x[i];
		}
		cons.add(v3<=nTestCap);
		v3.end();
		model.add(cons); 
	
		
	//---------------------Constraint Definition Ends here-------------------------------------------------------

	//---------------------Objective Function--------------------------------------------------------------------

		IloExpr obj(env);

		for(j=0;j<nFaults;j++){
			obj+=f[j];
		}

		for(j=0;j<nLines;j++){
			obj+=l[j];
		}
			
		model.add(IloMinimize(env, obj));			
			
	//---------------------Objective Function Def ends here----------------------------------------------------------
		
	//---------------------Solving the model-------------------------------------------------------------------------  
			
		IloCplex cplex(model);
		//if (easyFlag == 1) cplex.setParam(IloCplex::EpGap, 0.8); 
		//else 
		//cplex.setParam(IloCplex::EpGap, 0.25); 
		//cplex.setParam(IloCplex::EpGap, 0.005); 
		//if (easyFlag == 1) cplex.setParam(IloCplex::TiLim, 10);
		//else cplex.setParam(IloCplex::TiLim, 120);
		//else 
		cplex.setParam(IloCplex::TiLim, 21600);
		//cplex.setOut(env.getNullStream());
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);    //Select barrier algorithm 
        //cplex.setParam(IloCplex::BarCrossAlg, IloCplex::NoAlg);  //Without crossover

		//cplex.exportModel("lpex.lp");
		cplex.solve();
		
	//----------------------Extracting results for output-------------------------------------------------------------	
		
	tempDouble = cplex.getObjValue();
	
	fprintf(outfile, "Obj. Fun. %lf\n",cplex.getObjValue());  
	fprintf(outfile, "Gap. %lf\n", cplex.getMIPRelativeGap());

	////fprintf(outfile, "%lf ",cplex.getObjValue());  

	fprintf(outfile, "Selected Tests\n"); 
	
	nSelectedTests = 0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		//if (temp>0) fprintf (outfile,"%d %s\n", j+1, test_list[j].code);  
		nSelectedTests += temp; 
	}
	
	fprintf(outfile, "Number of Selected Tests %d\n", nSelectedTests); 
	fprintf(outfile2, "%d ", nSelectedTests); 
	
	memset(selectedTests, 0, (unsigned)(nSelectedTests) * (unsigned)sizeof(int));
	temp2=0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		if (temp>0) {
			selectedTests[temp2] = j; 
			temp2 = temp2 + 1;  
		}
	}


	fprintf(outfile, "Checking Line Coverage\n"); 

	temp2=0; 
	for(j=0;j<nLines;j++){
		if (lineFlag[j]>0){
			temp = 0; 
			for(i=0;i<nSelectedTests;i++){
				//fprintf (outfile,"%d %s\n", selectedTests[i], test_list[selectedTests[i]].code);  
				for(k=0;k<test_list[selectedTests[i]].nLinecover;k++){
					if (test_list[selectedTests[i]].lines[k] ==j+1) temp =1; 
				}
			}
			if (temp <1) fprintf(outfile, "Line Not Covered %d\n", j+1); 
			temp2 = temp2 + temp; 
		}
	}

	fprintf(outfile, "Number of Lines Covered %d\n", temp2); 
	fprintf(outfile2, "%d ", temp2); 
	LC = temp2; 
	objF = temp2; 

	fprintf(outfile, "Faults\n"); 
	
	temp2 = 0; 
	for(j=0;j<nFaults;j++){
		temp = cplex.getValue(f[j]);
		if (cplex.getValue(f[j])-temp > 0.8) temp = temp+1;
		if (temp>0) fprintf (outfile,"%d\n", faultID[j]);  
		temp2 += temp; 
	}
	
	fprintf(outfile, "Uncovered Faults %d\n", temp2); 
	fprintf(outfile2, "%d ", nFaults-temp2); 
	FC = nFaults-temp2; 
	objF += nFaults-temp2; 
	fprintf(outfile2, "%d ", objF); 

	for(j=0;j<nFaults;j++){
		temp = 0;
		for(i=0;i<nSelectedTests;i++){
			for(k=0;k<test_list[selectedTests[i]].nFaultcover;k++){
				if (test_list[selectedTests[i]].faults[k] == faultID[j]) temp = 1;
				//fprintf(outfile, "%d %d\n", test_list[selectedTests[i]].faults[k], temp); 
			}
		}
		if (temp <1) fprintf(outfile, "Fault Not Covered %d\n", faultID[j]); 
	}
		


	}

   catch (IloException& e) {
      cerr << "Error: " << e.getMessage() << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

	env.end();

	return tempDouble;
	


}

double Solver2(){

	int i,j,k,l,m;
	int temp, temp2, temp3; 
	double tempDouble; 
	int * selectedTests; 
	int nSelectedTests; 

		
	IloEnv env; //Intialize the Cplex Environment 
	try{
		IloModel model(env);
		
	//---------------------Decision variable Definition Starts here---------------------------------------------

		//IloNumVarArray Inv(env, planHorizon, 0, IloInfinity, ILOINT);
		//IloNumVarArray Prod(env, planHorizon, 0, IloInfinity, ILOINT);
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);

		//D2Array Inv(env,nItems); 
		//for(j=0;j<nItems;j++){
		//	Inv[j] = IloNumVarArray(env,planHorizon, 0, IloInfinity, ILOINT); 	
		//}
	
	//---------------------Decision variable Definition Ends here-----------------------------------------------

	//---------------------Constraint Definition Starts here----------------------------------------------------	
		
		IloRangeArray cons(env);
		for(j=0;j<nLines;j++){
			if (lineFlag[j]>0){
				IloExpr v1(env); 
				for(i=0;i<nTests;i++){
					temp = 0;
					for(k=0;k<test_list[i].nLinecover;k++){
						if (test_list[i].lines[k] == j+1) temp = 1;
					}
					if (temp) v1+=x[i];
				}
					
				cons.add(v1>=1);
				v1.end();
			}
		}
		model.add(cons); 

		for(j=0;j<nFaults;j++){
			IloExpr v2(env); 
			for(i=0;i<nTests;i++){
				temp = 0;
				for(k=0;k<test_list[i].nFaultcover;k++){
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2+=x[i];
			}
			v2+=f[j];
			cons.add(v2>=1);
			v2.end();
		}
		model.add(cons); 


		for(j=0;j<nTopTen;j++){
			if (lineFlag[topTenLines[j]]>0){
				IloExpr v3(env); 
				for(i=0;i<nTests;i++){
					temp = 0;
					for(k=0;k<test_list[i].nLinecover;k++){
						if (test_list[i].lines[k] == topTenLines[j]+1) temp = 1;
					}
					if (temp) v3+=x[i];
				}
					
				cons.add(v3>=topTenCover[j]);
				v3.end();
			}
		}
		model.add(cons); 

		
	//---------------------Constraint Definition Ends here-------------------------------------------------------

	//---------------------Objective Function--------------------------------------------------------------------

		IloExpr obj(env);

		for(j=0;j<nFaults;j++){
			obj+=f[j];
		}

		for(j=0;j<nTests;j++){
			obj+=x[j];
		}
			

		model.add(IloMinimize(env, obj));			
			
	//---------------------Objective Function Def ends here----------------------------------------------------------
		
	//---------------------Solving the model-------------------------------------------------------------------------  
			
		IloCplex cplex(model);
		//if (easyFlag == 1) cplex.setParam(IloCplex::EpGap, 0.8); 
		//else 
		//cplex.setParam(IloCplex::EpGap, 0.25); 
		//cplex.setParam(IloCplex::EpGap, 0.005); 
		//if (easyFlag == 1) cplex.setParam(IloCplex::TiLim, 10);
		//else cplex.setParam(IloCplex::TiLim, 120);
		//else 
		cplex.setParam(IloCplex::TiLim, 21600);
		//cplex.setOut(env.getNullStream());
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);    //Select barrier algorithm 
        //cplex.setParam(IloCplex::BarCrossAlg, IloCplex::NoAlg);  //Without crossover

		//cplex.exportModel("lpex.lp");
		cplex.solve();
		
	//----------------------Extracting results for output-------------------------------------------------------------	
		
	tempDouble = cplex.getObjValue();
	
	fprintf(outfile, "Obj. Fun. %lf\n",cplex.getObjValue());  
	fprintf(outfile, "Gap. %lf\n", cplex.getMIPRelativeGap());

	////fprintf(outfile, "%lf ",cplex.getObjValue());  

	fprintf(outfile, "Selected Tests\n"); 
	
	nSelectedTests = 0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		//if (temp>0) fprintf (outfile,"%d %s\n", j+1, test_list[j].code);  
		nSelectedTests += temp; 
	}
	
	fprintf(outfile, "Number of Selected Tests %d\n", nSelectedTests); 
	fprintf(outfile2, "%d ", nSelectedTests); 
	
	memset(selectedTests, 0, (unsigned)(nSelectedTests)* (unsigned)sizeof(int));
	temp2=0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		if (temp>0) {
			selectedTests[temp2] = j; 
			temp2 = temp2 + 1;  
		}
	}


	fprintf(outfile, "Checking Line Coverage\n"); 

	for(j=0;j<nLines;j++){
		if (lineFlag[j]>0){
			temp = 0; 
			for(i=0;i<nSelectedTests;i++){
				//fprintf (outfile,"%d %s\n", selectedTests[i], test_list[selectedTests[i]].code);  
				for(k=0;k<test_list[selectedTests[i]].nLinecover;k++){
					if (test_list[selectedTests[i]].lines[k] ==j+1) temp =1; 
				}
			}
			if (temp <1) fprintf(outfile, "Line Not Covered %d\n", j+1); 
		}
	}


	fprintf(outfile, "Faults\n"); 
	
	temp2 = 0; 
	for(j=0;j<nFaults;j++){
		temp = cplex.getValue(f[j]);
		if (cplex.getValue(f[j])-temp > 0.8) temp = temp+1;
		if (temp>0) fprintf (outfile,"%d\n", faultID[j]);  
		temp2 += temp; 
	}
	
	fprintf(outfile, "Uncovered Faults %d\n", temp2); 
	fprintf(outfile2, "%d ", nFaults-temp2); 
	fprintf(outfile2, "%d ", nSelectedTests+temp2); 

	for(j=0;j<nFaults;j++){
		temp = 0;
		for(i=0;i<nSelectedTests;i++){
			for(k=0;k<test_list[selectedTests[i]].nFaultcover;k++){
				if (test_list[selectedTests[i]].faults[k] == faultID[j]) temp = 1;
				//fprintf(outfile, "%d %d\n", test_list[selectedTests[i]].faults[k], temp); 
			}
		}
		if (temp <1) fprintf(outfile, "Fault Not Covered %d\n", faultID[j]); 
	}
		


	}

   catch (IloException& e) {
      cerr << "Error: " << e.getMessage() << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

	env.end();

	return tempDouble;
	


}

double Solver(){

	int i,j,k,l,m;
	int temp, temp2, temp3; 
	double tempDouble; 
	int * selectedTests; 
	int nSelectedTests; 

		
	IloEnv env; //Intialize the Cplex Environment 
	try{
		IloModel model(env);
		
	//---------------------Decision variable Definition Starts here---------------------------------------------

		//IloNumVarArray Inv(env, planHorizon, 0, IloInfinity, ILOINT);
		//IloNumVarArray Prod(env, planHorizon, 0, IloInfinity, ILOINT);
		IloNumVarArray x(env, nTests, 0, 1, ILOINT);
		IloNumVarArray f(env, nFaults, 0, 1, ILOINT);

		//D2Array Inv(env,nItems); 
		//for(j=0;j<nItems;j++){
		//	Inv[j] = IloNumVarArray(env,planHorizon, 0, IloInfinity, ILOINT); 	
		//}
	
	//---------------------Decision variable Definition Ends here-----------------------------------------------

	//---------------------Constraint Definition Starts here----------------------------------------------------	
		
		IloRangeArray cons(env);
		for(j=0;j<nLines;j++){
			if (lineFlag[j]>0){
				IloExpr v1(env); 
				for(i=0;i<nTests;i++){
					temp = 0;
					for(k=0;k<test_list[i].nLinecover;k++){
						if (test_list[i].lines[k] == j+1) temp = 1;
					}
					if (temp) v1+=x[i];
				}
					
				cons.add(v1>=1);
				v1.end();
			}
		}
		model.add(cons); 

		for(j=0;j<nFaults;j++){
			IloExpr v2(env); 
			for(i=0;i<nTests;i++){
				temp = 0;
				for(k=0;k<test_list[i].nFaultcover;k++){
					if (test_list[i].faults[k] == faultID[j]) temp = 1;
				}
				if (temp) v2+=x[i];
			}
			v2+=f[j];
			cons.add(v2>=1);
			v2.end();
		}
		model.add(cons); 


	
		
	//---------------------Constraint Definition Ends here-------------------------------------------------------

	//---------------------Objective Function--------------------------------------------------------------------

		IloExpr obj(env);

		for(j=0;j<nFaults;j++){
			obj+=f[j];
		}

		for(j=0;j<nTests;j++){
			obj+=x[j];
		}
			

		model.add(IloMinimize(env, obj));			
			
	//---------------------Objective Function Def ends here----------------------------------------------------------
		
	//---------------------Solving the model-------------------------------------------------------------------------  
			
		IloCplex cplex(model);
		//if (easyFlag == 1) cplex.setParam(IloCplex::EpGap, 0.8); 
		//else 
		//cplex.setParam(IloCplex::EpGap, 0.25); 
		//cplex.setParam(IloCplex::EpGap, 0.005); 
		//if (easyFlag == 1) cplex.setParam(IloCplex::TiLim, 10);
		//else cplex.setParam(IloCplex::TiLim, 120);
		//else 
		cplex.setParam(IloCplex::TiLim, 21600);
		//cplex.setOut(env.getNullStream());
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier);    //Select barrier algorithm 
        //cplex.setParam(IloCplex::BarCrossAlg, IloCplex::NoAlg);  //Without crossover

		//cplex.exportModel("lpex.lp");
		cplex.solve();
		
	//----------------------Extracting results for output-------------------------------------------------------------	
		
	tempDouble = cplex.getObjValue();
	
	fprintf(outfile, "Obj. Fun. %lf\n",cplex.getObjValue());  
	fprintf(outfile, "Gap. %lf\n", cplex.getMIPRelativeGap());

	////fprintf(outfile, "%lf ",cplex.getObjValue());  

	fprintf(outfile, "Selected Tests\n"); 
	
	nSelectedTests = 0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		//if (temp>0) fprintf (outfile,"%d %s\n", j+1, test_list[j].code);  
		nSelectedTests += temp; 
	}
	
	fprintf(outfile, "Number of Selected Tests %d\n", nSelectedTests); 
	fprintf(outfile2, "%d ", nSelectedTests); 

	memset(selectedTests, 0, (unsigned)(nSelectedTests) * (unsigned)sizeof(int));
	temp2=0; 
	for(j=0;j<nTests;j++){
		temp = cplex.getValue(x[j]);
		if (cplex.getValue(x[j])-temp > 0.8) temp = temp+1;
		if (temp>0) {
			selectedTests[temp2] = j; 
			temp2 = temp2 + 1;  
		}
	}


	fprintf(outfile, "Checking Line Coverage\n"); 

	for(j=0;j<nLines;j++){
		if (lineFlag[j]>0){
			temp = 0; 
			for(i=0;i<nSelectedTests;i++){
				//fprintf (outfile,"%d %s\n", selectedTests[i], test_list[selectedTests[i]].code);  
				for(k=0;k<test_list[selectedTests[i]].nLinecover;k++){
					if (test_list[selectedTests[i]].lines[k] ==j+1) temp =1; 
				}
			}
			if (temp <1) fprintf(outfile, "Line Not Covered %d\n", j+1); 
		}
	}


	fprintf(outfile, "Faults\n"); 
	
	temp2 = 0; 
	for(j=0;j<nFaults;j++){
		temp = cplex.getValue(f[j]);
		if (cplex.getValue(f[j])-temp > 0.8) temp = temp+1;
		if (temp>0) fprintf (outfile,"%d %d\n", j, faultID[j]);  
		temp2 += temp; 
	}
	
	fprintf(outfile, "Uncovered Faults %d\n", temp2); 
	fprintf(outfile2, "%d ", nFaults-temp2); 
	fprintf(outfile2, "%d ", nSelectedTests+temp2); 


	for(j=0;j<nFaults;j++){
		temp = 0;
		for(i=0;i<nSelectedTests;i++){
			for(k=0;k<test_list[selectedTests[i]].nFaultcover;k++){
				if (test_list[selectedTests[i]].faults[k] == faultID[j]) temp = 1;
				//fprintf(outfile, "%d %d\n", test_list[selectedTests[i]].faults[k], temp); 
			}
		}
		if (temp <1) fprintf(outfile, "Fault Not Covered %d\n", faultID[j]); 
	}
		


	}

   catch (IloException& e) {
      cerr << "Error: " << e.getMessage() << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

	env.end();

	return tempDouble;
	


}

char *trim(char *str)
{
    size_t len = 0;
    char *frontp = str;
    char *endp = NULL;

    if( str == NULL ) { return NULL; }
    if( str[0] == '\0' ) { return str; }

    len = strlen(str);
    endp = str + len;

    /* Move the front and back pointers to address the first non-whitespace
     * characters from each end.
     */
    while( isspace((unsigned char) *frontp) ) { ++frontp; }
    if( endp != frontp )
    {
        while( isspace((unsigned char) *(--endp)) && endp != frontp ) {}
    }

    if( str + len - 1 != endp )
            *(endp + 1) = '\0';
    else if( frontp != str &&  endp == frontp )
            *str = '\0';

    /* Shift the string so that it starts at str so that if it's dynamically
     * allocated, we can still free it on the returned pointer.  Note the reuse
     * of endp to mean the front of the string buffer now.
     */
    endp = str;
    if( frontp != str )
    {
            while( *frontp ) { *endp++ = *frontp++; }
            *endp = '\0';
    }


    return str;
}

int compare_ints (const void *a, const void *b){
	 if ((*(int **)a)[1] < (*(int **)b)[1]) return 1;
	 else return -1;
}

void readdata(){

	
	int i, j, k;
	int tempInt, tempInt1, tempInt2, tempInt3, tempInt4; 
	int lines; 
	char tempChar; 
	char s[200];
	int count = 0;
	int foundLetter = 0;
	int ch2=0; 
	int ch=0; 
	int * templist; 
	const char EOL = '\n';
	char str[15]; 
	int tempFault; 
	int faultFl; 
	int nLinesCovered; 
	double tempDouble; 
	
	nTests = 0; 
	while ((ch = fgetc(infile3)) != EOF){
		if (ch == '\n')
		nTests++;
    }
	fclose(infile3);

	fprintf (outfile,"Tests: %d\n",nTests);

	test_list = (struct test *) malloc(nTests * sizeof(struct test));

	infile3 = fopen("rtime.info","r");
	
	for(i=0; i<nTests ;i++) {
		test_list[i].id = i; 
		fscanf(infile3, "%s", &test_list[i].code);
		//fprintf(outfile, "%s ", test_list[i].code);
		fscanf(infile3, "%d", &test_list[i].time);
		//fprintf(outfile, "%d\n", test_list[i].time);
	}

	fclose(infile3); 

	tempInt = 0; 
	while ((ch = fgetc(infile)) != EOF){
		if (ch == '\n')
		tempInt++;
    }
	fclose(infile);

	fprintf (outfile,"Tests: %d\n",tempInt);
	if (tempInt!=nTests) fprintf (outfile,"DATA DO NOT MATCH\n");

	
	infile = fopen("cov.info","r");
	for(i=0; i<nTests ;i++) {
	//for(i=0; i<5 ;i++) {
		fscanf(infile, "%s", &str);
		tempInt = -1; 
		for(j=0; j<nTests ;j++) {
			if(strcmp(test_list[j].code,str) == 0) tempInt=j; 
		}
		
		//fprintf(outfile, "%d\n", tempInt);
		tempInt2 = 0; 
		//while (((ch = fgetc(infile)) != EOL) || ((ch = fgetc(infile)) != EOF)){
		while ((ch = fgetc(infile)) != EOL){

			if (ch == ' ')
			tempInt2++;
		}
		//fprintf (outfile,"Lines: %d\n",tempInt2);

		test_list[tempInt].nLinecover = tempInt2;

		test_list[tempInt].lines = (int *) malloc(test_list[tempInt].nLinecover * sizeof(int));

	}

	fclose(infile);


	
	infile = fopen("cov.info","r");
	for(i=0; i<nTests ;i++) {
	//for(i=0; i<5 ;i++) {
		fscanf(infile, "%s", &str);
		tempInt = -1; 
		for(j=0; j<nTests ;j++) {
			if(strcmp(test_list[j].code,str) == 0) tempInt=j; 
		}

		//fprintf(outfile, "%d\n", tempInt);

		for(j=0; j<test_list[tempInt].nLinecover ;j++) {
			fscanf(infile, "%d", &test_list[tempInt].lines[j]);
			//fprintf(outfile, "%d ", test_list[tempInt].lines[j]);
		}
		//fprintf (outfile,"\n");
	}

	fclose(infile);

	memset(lineFlag, 0, (unsigned)(nLines) * (unsigned)sizeof(int));
	for(i=0; i<nTests ;i++) {
		for(j=0; j<test_list[i].nLinecover ;j++) {
			//if(test_list[i].lines[j]<1) fprintf (outfile,"ALOO %d\n", i+1);
			//if(test_list[i].lines[j]>nLines) fprintf (outfile,"ALOO %d %d\n", i+1, test_list[i].lines[j]);
			lineFlag[test_list[i].lines[j]-1] = 1; 
		}
	}
	

	
	nLinesCovered = 0; 
	for(i=0; i<nLines ;i++) {
		if (lineFlag[i] <1) fprintf (outfile,"Line not covered %d\n", i+1);
		else nLinesCovered++;  
	}


	fprintf (outfile,"nLinesCovered %d\n", nLinesCovered);

	tempInt = 0; 
	while ((ch = fgetc(infile2)) != EOF){
		if (ch == '\n')
		tempInt++;
    }
	fclose(infile2);

	fprintf (outfile,"Tests: %d\n",tempInt);
	if (tempInt!=nTests) fprintf (outfile,"DATA DO NOT MATCH\n");

	memset(tempList, 0, (unsigned)(nTests)* (unsigned)sizeof(int));

	infile2 = fopen("fault.info","r");
	for(i=0; i<nTests ;i++) {
	//for(i=0; i<1500 ;i++) {
		//fscanf(infile2, "%s", &str);
		//tempInt = -1; 
		//for(j=0; j<nTests ;j++) {
		//	if(strcmp(test_list[j].code,str) == 0) tempInt=j; 
		//}

		//tempInt2 = 0; 
		//ch2=0;
		//while ((ch = fgetc(infile2)) != EOL){

		//	if ((ch == ' ') && (ch2 != ' '))
		//	tempInt2++;

		//	ch2 = ch; 
		//}
		//fprintf (outfile,"Lines: %d\n",tempInt2);

		char buffer[20000];
		int numbers[10000];
		char *ptr;
		int cnt = 0;

		//fprintf(outfile, "%d\n", tempInt);
		fgets(&buffer[0], 20000, infile2); // Get the string
		printf("%s\n", trim(buffer));
		ptr = strtok(buffer, " "); // Split in every space
		do {
            //fprintf(outfile, "Number %d: %s\n", cnt, ptr);
            numbers[cnt] = strtol(ptr, NULL, 0); // Grab number
            if(strcmp(ptr," ") == 0)  fprintf(outfile, "Space\n");
			else cnt++;
		} while((ptr = strtok(NULL, " ")));
		//fprintf(outfile, "Total numbers: %d\n", cnt-1);
		templist[i] = cnt-1; 

		//const int LENGTH_LINE = 100;
		//char line[LENGTH_LINE];
		//int len;
	

		//fgets(line,LENGTH_LINE,infile3);
		//len = strlen(line);
		//if(line[len-1] == '\n')
		//	fprintf(outfile,"I've a line\n");
		//	
		//fprintf(outfile,"%s\n", line);
		//fprintf (outfile,"Faults: %d\n",tempInt2);

		//test_list[tempInt].nFaultcover = tempInt2;

		//test_list[tempInt].faults = (int *) malloc(test_list[tempInt].nFaultcover * sizeof(int));

	}

	fclose(infile2);
		
	memset(faultID, 0, (unsigned)(nFaults)* (unsigned)sizeof(int));
	infile2 = fopen("fault.info","r");
	tempFault = 0; 
	for(i=0; i<nTests ;i++) {
	//for(i=0; i<5 ;i++) {
		fscanf(infile2, "%s", &str);
		tempInt = -1; 
		for(j=0; j<nTests ;j++) {
			if(strcmp(test_list[j].code,str) == 0) tempInt=j; 
		}

		//fprintf(outfile, "%d\n", tempInt);

		test_list[tempInt].nFaultcover = templist[i]; 

		test_list[tempInt].faults = (int *) malloc(test_list[tempInt].nFaultcover * sizeof(int));

		for(j=0; j<test_list[tempInt].nFaultcover ;j++) {
			fscanf(infile2, "%d", &test_list[tempInt].faults[j]);
			faultFl = 0; 
			for(k=0; k<tempFault ;k++) {
				if(faultID[k]==test_list[tempInt].faults[j]) faultFl = 1;
			}
			if (faultFl<1) {
				faultID[tempFault] = test_list[tempInt].faults[j]; 
				tempFault += 1; 
			}
			///fprintf(outfile, "%d ", test_list[tempInt].faults[j]);
		}
		//fprintf (outfile,"\n");
	}

	fclose(infile);

	fprintf (outfile,"tempFault %d\n", tempFault);

	nTopTen = nLinesCovered*0.1;
	//tempDouble = (double) (nLinesCovered*0.1);
	//if (tempDouble -nTopTen>0.5) nTopTen +=1;  

	fprintf (outfile,"TopTen %d\n", nTopTen);

	memset(topTenLines, 0, (unsigned)(nTopTen)* (unsigned)sizeof(int));
	memset(topTenCover, 0, (unsigned)(nTopTen)* (unsigned)sizeof(int));
	memset(lineCover, 0, (unsigned)(nLines)* (unsigned)sizeof(int));
	for(i=0; i < nLines; i++){
		memset(lineCover[i], 0, 2 * (unsigned)sizeof(int));
	}

	for(i=0; i<nTests ;i++) {
		for(j=0; j<test_list[i].nLinecover ;j++) {
			lineCover[test_list[i].lines[j]-1][1] += 1; 
		}
	}

	for(i=0; i<nLines ;i++) {
		lineCover[i][0] = i; 
		//fprintf (outfile,"%d\n", lineCover[i][1]);
	}

	qsort(lineCover, nLines, sizeof(int*), compare_ints);

	for(i=0; i<nTopTen ;i++) {
		//fprintf (outfile,"%d %d\n", lineCover[i][0], lineCover[i][1]);
		topTenLines[i] = lineCover[i][0];
		topTenCover[i] = 0.1*lineCover[i][1];
		//tempDouble = (double) (0.1*lineCover[i][1]);
		//if (tempDouble - topTenCover[i]>0.5) topTenCover[i] +=1;  
		//fprintf (outfile,"%d %d\n", topTenLines[i], topTenCover[i]);
	}

}





void main(int argc, char *argv[]){

	int i,j; 
	int Budget; 
	double tempDouble; 

	time_t seconds;
	time_t seconds2;
	time_t seconds3;
	time_t seconds4;
	time_t seconds5;
	time_t seconds6;
	time_t seconds7;
	time_t seconds8;

	nLines = atoi(argv[1]); ; 
	nFaults = atoi(argv[2]); ; 
	//nLines = 3143; 
	//nLines = 26051;
	//nFaults = 41631; 


	infile = fopen("cov.info","r");
	infile2 = fopen("fault.info","r");
	infile3 = fopen("rtime.info","r");
	//infile4 = fopen("Classroom.dat","r");
	outfile = fopen("out.dat","w");
	outfile2 = fopen("sum.dat","w");

	readdata(); 


	seconds = time (NULL);
	Solver();
	seconds2 = time (NULL);
	fprintf (outfile, "Model 1 Solution Time (secs): %ld\n", (seconds2-seconds));
	fprintf (outfile2, "%ld\n", (seconds2-seconds));
	fprintf (outfile,"\n\n");

	
	Solver2();
	seconds3 = time (NULL);
	fprintf (outfile, "Model 2 Solution Time (secs): %ld\n", (seconds3-seconds2));
	fprintf (outfile2, "%ld\n", (seconds3-seconds2));
	fprintf (outfile,"\n\n");

	Budget = nTests*0.05;
	tempDouble = (double) (nTests*0.05);
	if (tempDouble -Budget>0.5) Budget +=1;  
	Solver3(Budget);
	seconds4 = time (NULL);
	fprintf (outfile, "Model 3-1 Solution Time (secs): %ld\n", (seconds4-seconds3));
	fprintf (outfile2, "%ld\n", (seconds4-seconds3));
	fprintf (outfile,"\n\n");
	
	Budget = nTests*0.1;
	tempDouble = (double) (nTests*0.1);
	if (tempDouble -Budget>0.5) Budget +=1;  
	Solver3(Budget);
	seconds5 = time (NULL);
	fprintf (outfile, "Model 3-2 Solution Time (secs): %ld\n", (seconds5-seconds4));
	fprintf (outfile2, "%ld\n", (seconds5-seconds4));
	fprintf (outfile,"\n\n");

	Budget = nTests*0.15;
	tempDouble = (double) (nTests*0.15);
	if (tempDouble -Budget>0.5) Budget +=1;  
	Solver3(Budget);
	seconds6 = time (NULL);
	fprintf (outfile, "Model 3-3 Solution Time (secs): %ld\n", (seconds6-seconds5));
	fprintf (outfile2, "%ld\n", (seconds6-seconds5));
	fprintf (outfile,"\n\n");

	Budget = nTests*0.2;
	tempDouble = (double) (nTests*0.2);
	if (tempDouble -Budget>0.5) Budget +=1;  
	Solver3(Budget);
	seconds7 = time (NULL);
	fprintf (outfile, "Model 3-4 Solution Time (secs): %ld\n", (seconds7-seconds6));
	fprintf (outfile2, "%ld\n", (seconds7-seconds6));

	Solver4(Budget);
	seconds8 = time (NULL);
	fprintf (outfile, "Model 4 Solution Time (secs): %ld\n", (seconds8-seconds7));
	fprintf (outfile2, "%ld\n", (seconds8-seconds7));
	

	fclose(outfile); 
	
}


