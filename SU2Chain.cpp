#include <iostream>
#include "config.h"
#include <stdlib.h>
#include<unistd.h>
#include <string>
#include <fstream>
#include <streambuf>
#include<armadillo>
#include<cassert>
#include<filesystem>

using namespace arma;


class SU2Driver {


public:


	rowvec designVariables;
	int numberOfDesignVariables;

	std::string filenameDesign = "dv.dat";
	std::string filenameForCL = "outputCL.dat";
	std::string filenameForCD = "outputCD.dat";
	std::string filenameForArea = "outputArea.dat";
	std::string filenameForIterationNumber = "iterationNumber.dat";
	std::string filenameGEOInput = "of_eval.csv";
	std::string filenameHistory  = "history_direct.csv";
	std::string filenameCDGradients  = "of_grad_cd.csv";
	std::string filenameCLGradients  = "of_grad_cl.csv";
	std::string filenameDoEInput = "DoEInput.csv";

	std::string workingDirectory;

	std::string directoryName;

	unsigned int iterationNumber = 0;

	mat dataInput;

	vec lb;
	vec ub;

	unsigned int numberOfDoESamples = 0;

	struct{

		double CL;
		double CD;
		vec gradCL;
		vec gradCD;

	} aerodynamicProperties;

	struct{

		double area;


	} geometricProperies;


	void changeToWorkingDirectory(void){

		assert(!workingDirectory.empty());
		chdir(workingDirectory.c_str());

	}


	void readDesignFromFile(void){

		assert(numberOfDesignVariables>0);

		designVariables = zeros<rowvec>(numberOfDesignVariables);
		std::fstream inputfile;
		inputfile.open(filenameDesign,ios::in);
		if (!inputfile)
			throw std::system_error(errno, std::system_category(), "failed to open "+filenameDesign);

		for(unsigned int i=0; i<numberOfDesignVariables; i++){
			inputfile>>designVariables(i);
		}

		designVariables.print("design variables = \n");
		inputfile.close();
	}

	void writeDesignToFile(void){

		std::fstream inputfile;
		inputfile.open(filenameDesign,ios::out);

		for(unsigned int i=0; i<numberOfDesignVariables; i++){
			inputfile<<designVariables(i)<<"\n";
		}
		inputfile.close();
	}


	void readIterationNumberFromFile(void){

		std::fstream inputfile;
		inputfile.open(filenameForIterationNumber,ios::in);
		inputfile>>iterationNumber;
		inputfile.close();


		iterationNumber++;


		std::cout<<"Iteration number = "<<iterationNumber<<"\n";

		inputfile.open(filenameForIterationNumber,ios::out);
		inputfile<<iterationNumber;
		inputfile.close();

	}

	void writeOutput(void){
		std::fstream inputfile;
		inputfile.open(filenameForCL,ios::out);
		inputfile<<aerodynamicProperties.CL;
		inputfile.close();
		inputfile.open(filenameForCD,ios::out);
		inputfile<<aerodynamicProperties.CD;
		inputfile.close();
		inputfile.open(filenameForArea,ios::out);
		inputfile<<geometricProperies.area;
		inputfile.close();

	}

	void writeOutputWithGradients(void){
		std::fstream inputfile;
		inputfile.open(filenameForCL,ios::out);
		inputfile<<aerodynamicProperties.CL<<"\n";

#if CL_GRADIENT_IS_ON
		for(unsigned int i=0; i<numberOfDesignVariables; i++){
			inputfile<<aerodynamicProperties.gradCL(i)<<"\n";
		}
#endif

		inputfile.close();
		inputfile.open(filenameForCD,ios::out);
		inputfile<<aerodynamicProperties.CD<<"\n";

		for(unsigned int i=0; i<numberOfDesignVariables; i++){
			inputfile<<aerodynamicProperties.gradCD(i)<<"\n";

		}

		inputfile.close();
		inputfile.open(filenameForArea,ios::out);
		inputfile<<geometricProperies.area;
		inputfile.close();

	}



	void clearFilesAndDirectories(void){

		system("rm DOE_* -rf");
		system("rm CL.csv");
		system("rm CD.csv");
		system("rm area.csv");

	}


	void createDirectoryDOE(unsigned int sampleID){

		assert(sampleID > 0);
		directoryName = "DOE_" + std::to_string(sampleID);
		std::filesystem::create_directory(directoryName);

	}

	void createDirectoryForOptimizationIteration(void){

		directoryName = "DSN_" + std::to_string(iterationNumber);
		std::filesystem::create_directory(directoryName);

	}


	void copyMeshAndConfigFiles(void){

		system("cp ../config_GEO.cfg ./");
		system("cp ../config_DEF_BASE.cfg ./");
		system("cp ../config_CFD.cfg ./");
		system("cp ../mesh_NACA0012_inv.su2 ./");
		system("cp ../naca0012geometry.csv ./");
	}

	void addSampleToData(std::string output){

		assert(numberOfDesignVariables > 0);


		if(output == "lift"){
			rowvec newSample(numberOfDesignVariables+1);
			for(unsigned int i=0; i<numberOfDesignVariables; i++){
				newSample(i) = designVariables(i);
			}
			newSample(numberOfDesignVariables) = aerodynamicProperties.CL;
			appendRowVectorToCSVData(newSample,"CL.csv");
		}

		if(output == "liftWithGradient"){
			rowvec newSample(2*numberOfDesignVariables+1);
			for(unsigned int i=0; i<numberOfDesignVariables; i++){
				newSample(i) = designVariables(i);
			}
			newSample(numberOfDesignVariables) = aerodynamicProperties.CL;
			for(unsigned int i=numberOfDesignVariables+1; i<2*numberOfDesignVariables+1; i++){
				newSample(i) = aerodynamicProperties.gradCL(i-numberOfDesignVariables-1);
			}
			appendRowVectorToCSVData(newSample,"CL.csv");
		}

		if(output == "drag"){
			rowvec newSample(numberOfDesignVariables+1);
			for(unsigned int i=0; i<numberOfDesignVariables; i++){
				newSample(i) = designVariables(i);
			}
			newSample(numberOfDesignVariables) = aerodynamicProperties.CD;
			appendRowVectorToCSVData(newSample,"CD.csv");
		}

		if(output == "dragWithGradient"){
			rowvec newSample(2*numberOfDesignVariables+1);
			for(unsigned int i=0; i<numberOfDesignVariables; i++){
				newSample(i) = designVariables(i);
			}
			newSample(numberOfDesignVariables) = aerodynamicProperties.CD;
			for(unsigned int i=numberOfDesignVariables+1; i<2*numberOfDesignVariables+1; i++){
				newSample(i) = aerodynamicProperties.gradCD(i-numberOfDesignVariables-1);
			}
			appendRowVectorToCSVData(newSample,"CD.csv");
		}



		if(output == "area"){
			rowvec newSample(numberOfDesignVariables+1);
			for(unsigned int i=0; i<numberOfDesignVariables; i++){
				newSample(i) = designVariables(i);
			}
			newSample(numberOfDesignVariables) = geometricProperies.area;
			appendRowVectorToCSVData(newSample,"area.csv");
		}



	}

	void appendRowVectorToCSVData(rowvec v, std::string filename){

		assert(v.size() > 0);
		std::ofstream outfile;

		outfile.open(filename, std::ios_base::app); // append instead of overwrite

		outfile.precision(10);
		for(unsigned int i=0; i<v.size()-1; i++){

			outfile << v(i) <<",";
		}
		outfile << v(v.size()-1)<<"\n";
		outfile.close();

	}




	void setUpperAndLowerBounds(double down, double up){

		assert(numberOfDesignVariables > 0);
		lb = zeros<vec>(numberOfDesignVariables);
		ub = zeros<vec>(numberOfDesignVariables);
		assert(up>down);
		lb.fill(down);
		ub.fill(up);
	}


	void generateDoEInputMatrix(void){

		assert(numberOfDoESamples>0);
		assert(numberOfDesignVariables > 0);
		assert(lb.size() > 0);
		assert(ub.size() > 0);

		dataInput = zeros<mat>(numberOfDoESamples, numberOfDesignVariables);

		for(unsigned int i=0; i<numberOfDoESamples; i++){

			rowvec randomVec(numberOfDesignVariables, fill::randu);

			for(unsigned int j=0; j<numberOfDesignVariables; j++){
				randomVec(j) = randomVec(j)*(ub(j)- lb(j)) + lb(j);
			}

			if(i==0) randomVec.fill(0.0);


			dataInput.row(i) = randomVec;

		}
#if 0
		dataInput.print();
#endif

	}


	void generateDoEInputMatrixAroundASample(vec sample, double epsilon){

		assert(numberOfDoESamples>0);
		assert(numberOfDesignVariables > 0);
		assert(lb.size() > 0);
		assert(ub.size() > 0);
		assert(epsilon > 0.0);

		vec deltaX = (ub - lb)*epsilon;
		vec lbAroundASample = sample - deltaX;
		vec ubAroundASample = sample + deltaX;

		for(unsigned int j=0; j<numberOfDesignVariables; j++){

			if(lbAroundASample(j) < lb(j)) lbAroundASample(j) = lb(j);
			if(ubAroundASample(j) > ub(j)) ubAroundASample(j) = ub(j);
		}

		dataInput = zeros<mat>(numberOfDoESamples, numberOfDesignVariables);

		for(unsigned int i=0; i<numberOfDoESamples; i++){

			rowvec randomVec(numberOfDesignVariables, fill::randu);

			for(unsigned int j=0; j<numberOfDesignVariables; j++){
				randomVec(j) = randomVec(j)*(ubAroundASample(j)- lbAroundASample(j)) + lbAroundASample(j);
			}

			dataInput.row(i) = randomVec;

		}
#if 0
		dataInput.print();
#endif

	}



	void readDoEInputMatrix(void){

		mat data;
		data.load( filenameDoEInput, csv_ascii );

		assert(data.n_cols == numberOfDesignVariables);

		dataInput = data;


	}


	void callGeo(void){

		system("SU2_GEO config_GEO.cfg > geo_output.dat");
	}

	void readGeoOutput(void){
		field<std::string> header;
		mat data;
		data.load( csv_name(filenameGEOInput, header) );

		geometricProperies.area = data(0,0);

		std::cout<<"area = "<< geometricProperies.area << "\n";

	}


	void readGradientVectors(void){

		field<std::string> header;
		mat data;
		bool ok = data.load( csv_name(filenameCDGradients, header) );

		if(ok == false)
		{
			cout << "Problem with loading " << filenameCDGradients << endl;
		}

		vec grad = data.col(1);
		aerodynamicProperties.gradCD = grad;

		trans(grad).print("gradient of CD");

		//		ok = data.load( csv_name(filenameCLGradients, header) );
		//
		//		if(ok == false)
		//		{
		//			cout << "Problem with loading " << filenameCLGradients << endl;
		//		}
		//
		//
		//		grad = data.col(1);
		//		aerodynamicProperties.gradCL = grad;


	}

	void callSU2PrimalSolver(void){
		system("SU2_CFD config_CFD.cfg > SU2_primal_output.dat");
	}

	void callSU2AdjointSolver(void){
		system("SU2_CFD_AD config_CFD_AD.cfg");
	}

	void callMeshDeformation(void){
		system("SU2_DEF config_DEF.cfg > meshdefo_output.dat");
	}

	void plotAirfoil(void){
		system("python ../../../plotAirfoil.py");
	}


	void readOutputSU2Primal(void){

		field<std::string> header;
		mat data;
		data.load( csv_name("history_direct.csv", header) );

		//		data.print();
		vec cl = data.col(9);
		aerodynamicProperties.CL = cl(cl.size()-1);
		vec cd = data.col(8);
		aerodynamicProperties.CD = cd(cd.size()-1);
#if 1
		std::cout<<"CD = "<<aerodynamicProperties.CD<<"\n";
		std::cout<<"CL = "<<aerodynamicProperties.CL<<"\n";
#endif
	}

	void callSU2SimulationChain(void){


		assert(designVariables.size() == numberOfDesignVariables);

		std::ifstream baseInput("config_DEF_BASE.cfg");
		std::string baseString((std::istreambuf_iterator<char>(baseInput)),
				std::istreambuf_iterator<char>());

		//		std::cout<<baseString;

		baseString +="DV_VALUE = ";
		for(int i=0;i<numberOfDesignVariables-1;i++){

			baseString +=std::to_string(designVariables(i));
			baseString +=",";
		}
		baseString +=std::to_string(designVariables(numberOfDesignVariables-1));
		baseString +="\n";
		//		std::cout<<baseString;
		std::ofstream out("config_DEF.cfg");
		out <<baseString;
		out.close();

		callMeshDeformation();
		callSU2PrimalSolver();
		readOutputSU2Primal();

		callGeo();
		readGeoOutput();


	}

	void callSU2AdjointChain(void){


		system("python script.py > adjoint_chain_output");
		system("cp ./DESIGNS/DSN_001/DIRECT/history_direct.csv ./");
		system("cp ./DESIGNS/DSN_001/ADJOINT_DRAG/of_grad_cd.csv ./");
		//		system("cp ./DESIGNS/DSN_001/ADJOINT_LIFT/of_grad_cl.csv ./");
		system("cp ./DESIGNS/DSN_001/GEOMETRY/of_eval.csv ./");
		readOutputSU2Primal();
		readGeoOutput();
		readGradientVectors();
		writeOutputWithGradients();

	}



	void performOptimizationIterationWithoutAdjoints(void){

		std::cout<<"Calling simulation chain...\n";
		readIterationNumberFromFile();
		readDesignFromFile();
		createDirectoryForOptimizationIteration();
		std::string path = "./" + directoryName;
		//			std::cout<<"path = "<<path<<"\n";
		chdir(path.c_str());
		copyMeshAndConfigFiles();
		std::cout<<"Calling SU2_CFD...\n";
		callSU2SimulationChain();
		plotAirfoil();
		path = "../";
		chdir(path.c_str());
		writeOutput();

	}

	void performOptimizationIterationWithAdjoints(void){

		std::cout<<"Calling simulation chain with adjoints ...\n";

		std::string path = "./SimulationData";
		chdir(path.c_str());
		std::string  command = "cp ../dv.dat ./";
		system(command.c_str());

		readIterationNumberFromFile();
		createDirectoryForOptimizationIteration();
		callSU2AdjointChain();

		command = "cp ./DESIGNS/DSN_001/* ./";
		command+= directoryName;
		command+= " -r";
//		std::cout<<command<<std::endl;
		system(command.c_str());
		command = "cp outputCD.dat ../";
		system(command.c_str());
		command = "cp outputCL.dat ../";
		system(command.c_str());
		command = "cp outputArea.dat ../";
		system(command.c_str());
		path = "../";
		chdir(path.c_str());

	}



	vec evaluateGradientWithFiniteDifferences(void){

		assert(numberOfDesignVariables>0);
		assert(designVariables.size() > 0);

		std::cout<<"Evaluating finite differences for the baseline = \n";
		designVariables.print();
		vec grad(numberOfDesignVariables);

		rowvec dvSave = designVariables;

		for(unsigned int i=0; i<numberOfDesignVariables; i++){

			std::cout<<"Finite difference for the variable number = "<<i<<"\n";

			double epsilon = designVariables(i)*0.0000001;

			if(epsilon < 10E-14){

				epsilon = 10E-6;
			}

			designVariables(i) += epsilon;

			callSU2SimulationChain();

			double fplus = aerodynamicProperties.CD;

			std::cout<<"fplus = "<<fplus<<"\n";


			designVariables(i) -= 2*epsilon;
			callSU2SimulationChain();

			double fminus = aerodynamicProperties.CD;

			std::cout<<"fminus = "<<fminus<<"\n";

			double fd = (fplus - fminus)/(2*epsilon);
			grad(i) = fd;
			std::cout<<"finite difference = "<<fd<<"\n";

			designVariables = dvSave;


		}


		return grad;
	}



	void performDoE(void){

		assert(numberOfDoESamples > 0);
		assert(dataInput.n_rows > 0);

		clearFilesAndDirectories();

		for(unsigned int i=0; i<numberOfDoESamples; i++){

			std::cout<<"DOE Sample #"<<i<<"\n";

			unsigned int sampleID = i+1;
			createDirectoryDOE(sampleID);
			workingDirectory = "./" + directoryName;
			changeToWorkingDirectory();
			copyMeshAndConfigFiles();
			designVariables = dataInput.row(i);
			callSU2SimulationChain();
			std::cout<<"Simulation is complete...\n";
			plotAirfoil();
			workingDirectory = "../";
			changeToWorkingDirectory();
			addSampleToData("lift");
			addSampleToData("drag");
			addSampleToData("area");


		}





	}

	void performDoEWithAdjoints(void){

		clearFilesAndDirectories();

		//		readDoEInputMatrix();
		//		dataInput.print();
		numberOfDoESamples = dataInput.n_rows;

		for(unsigned int i=0; i<numberOfDoESamples; i++){

			std::cout<<"DOE Sample #"<<i<<"\n";

			unsigned int sampleID = i+1;
			createDirectoryDOE(sampleID);


			designVariables = dataInput.row(i);
			designVariables.print("design variables =\n");
			writeDesignToFile();
			callSU2AdjointChain();

			addSampleToData("lift");
			addSampleToData("dragWithGradient");
			addSampleToData("area");

			std::string command = "cp ./DESIGNS/DSN_001/* ./";
			command+= directoryName;
			command+= " -r";
			std::cout<<command<<std::endl;
			system(command.c_str());

		}


	}


};



int main(int argc, char **argv) {

	SU2Driver driver;
	unsigned int dim = 20;
	driver.numberOfDesignVariables = dim;
	driver.setUpperAndLowerBounds(-0.001,0.001);

	if(argc == 1){

		std::cout<<"ERROR: working mode is not specified\n";
		std::cout<<"Example call:\n";
		std::cout<<"./SU2Chain primal_design\n";
		abort();
	}

	std::string mode = argv[1];
	std::cout << "Working mode = " << mode<<std::endl;



	if(mode == "primal_design"){


		driver.workingDirectory = "./naca0012/OPTIMIZATION_WO_ADJOINTS";
		driver.changeToWorkingDirectory();

		driver.performOptimizationIterationWithoutAdjoints();

	}

	if(mode == "adjoint_design"){


		driver.workingDirectory = "./naca0012/OPTIMIZATION_WITH_ADJOINTS";
		driver.changeToWorkingDirectory();
		driver.performOptimizationIterationWithAdjoints();

	}

	if(mode == "primal_DoE"){


		vec dv(driver.numberOfDesignVariables, fill::zeros);
		driver.numberOfDoESamples = 100;
		driver.generateDoEInputMatrixAroundASample(dv,1.0);
		driver.workingDirectory = "./naca0012/DOE";
		driver.changeToWorkingDirectory();
		driver.performDoE();

	}


	if(mode == "adjoint_DoE"){


		vec dv(driver.numberOfDesignVariables, fill::zeros);
		driver.numberOfDoESamples = 100;
		driver.generateDoEInputMatrix();
		driver.workingDirectory = "./naca0012/DOE_ADJOINT";
		driver.changeToWorkingDirectory();
		driver.performDoEWithAdjoints();

	}


	return 0;
}
