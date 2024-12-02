//simplified grid class 
//member data GridID, Task, Rank, Dimension,Start/End Index, Left/Right Edge, Time 
//StaticSubgrids, NumberOfBaryonFields, FieldType, BaryonFileName, CourantSafetyNumber
//PPM params, NextGridThisLevel pointer, NextGridNextLevel pointer
//
//All functions were tested with mpi implementation w/ 2 Nodes and 4 processes
//.cpp files for parse_hierarchy, distribute_grids, etc. are test programs for functions 
//included here
//
//May break without mpi 
//subject to change, specifically opemmp implementation
//

#include<iostream>
#include<cstring>
using namespace std;

class grid
{
	private: //member data 
		int GridID; 
		int Task;
		int GridRank; 
		int GridDimension[3]; 
		int GridStartIndex[3];
		int GridEndIndex[3]; 
		float GridLeftEdge[3]; 
		float GridRightEdge[3]; 
		int GridLevel; 
		float Time; 
		int SubgridsAreStatic; 
		int NumberOfBaryonFields; 
		int FieldType[100]; 
		char BaryonFileName[100]; 
		float CourantSafetyNumber; 
		int PPMFlatteningParameter; 
		int PPMDiffusionParameter; 
		int PPMSteepeningParameter; 
		int NumberOfParticles; 
		int NumberOfActiveParticles; 
		int GravityBoundaryType; 
		int *NextGridThisLevel; 
		int *NextGridNextLevel;
		float *BaryonField[20]; //pointers to BaryonFields 


	public: //member functions
		grid(); //constructor which sets everything to defaults

		//setter routines 
		int SetGridID(int myGridID);
		int SetTask(int myTask);
		int SetGridRank(int myGridRank);
		int SetGridDimension(int myGridDimension[3]); 
		int SetGridStartIndex(int myGridStartIndex[3]);
		int SetGridEndIndex(int myGridEndIndex[3]); 
		int SetGridLeftEdge(float myGridLeftEdge[3]); 
		int SetGridRightEdge(float myGridRightEdge[3]);
		int SetGridLevel(int myGridLevel); 
		int SetTime(float myTime); 
		int SetSubgridsAreStatic(int mySubgridsAreStatic); 
		int SetNumberOfBaryonFields(int myNumberOfBaryonFields); 
		int SetFieldType(int myFieldType[100]); 
		int SetBaryonFileName(char *myBaryonFileName); 
		int SetCourantSafetyNumber(float myCourantSafetyNumber); 
		int SetPPMFlatteningParameter(int myPPMFlatteningParameter); 
		int SetPPMDiffusionParameter(int myPPMDiffusionParameter); 
		int SetPPMSteepeningParameter(int myPPMSteepeningParameter); 
		int SetNumberOfParticles(int myNumberOfParticles); 
		int SetNumberOfActiveParticles(int myNumberOfActiveParticles);
	    int SetGravityBoundaryType(int myGravityBoundaryType); 
		int SetNextGridThisLevel(int *myNextGridThisLevel); 
		int SetNextGridNextLevel(int *myNextGridNextLevel);

		//Accessor Routines
		int GetGridID();
	    int GetTask();	
		int GetGridRank(); 
		int* GetGridDimension();	
		int* GetGridStartIndex(); 
		int* GetGridEndIndex(); 
		float* GetGridLeftEdge(); 
		float* GetGridRightEdge(); 
		int GetGridLevel(); 
		float GetTime(); 
		int GetSubgridsAreStatic(); 
		int GetNumberOfBaryonFields(); 
		int* GetFieldType(); 
		char* GetBaryonFileName(); 
		float GetCourantSafetyNumber(); 
		int GetPPMFlatteningParameter(); 
		int GetPPMDiffusionParameter(); 
		int GetPPMSteepeningParameter(); 
		int GetNumberOfParticles(); 
		int GetNumberOfActiveParticles(); 
		int GetGravityBoundaryType();
		int* GetNextGridThisLevel(); 
		int* GetNextGridNextLevel();
	    void PrintAllData();	
		
		//friend functions
		friend int parse_hierarchy(const char *filename, int lower_grid, int upper_grid, grid localgrid[],int num_grids, int nodenum); 
		friend int ReadListOfInts(FILE *fptr, int N, int num[], int TestGridID);
		friend int distribute_grids(grid grids[],int num_nodes, int global_rank, const char *filename);
};

grid::grid()
{
	//set values to either zero or null, whichever makes sense
	GridID = 0; 
	Task = 0;
	GridRank = 0;
        GridLevel = 0; 
	Time = 0.0; 
	SubgridsAreStatic = 0; 
	NumberOfBaryonFields = 0; 
	CourantSafetyNumber = 0.0; 
	PPMFlatteningParameter = 0; 
	PPMDiffusionParameter = 0; 
	PPMSteepeningParameter = 0; 
	NumberOfParticles = 0; 
	NumberOfActiveParticles = 0; 
	GravityBoundaryType = 0; 
	for(int i = 0; i < 3; i++){
		GridDimension[i] = 0; 
		GridStartIndex[i] = 0; 
		GridEndIndex[i] = 0; 
		GridLeftEdge[i] = 0.0; 
		GridRightEdge[i] = 0.0;
	}
	for(int i = 0; i < 100; i++){
		FieldType[i] = 104; //This is what field undefined is in enzo's grid so following that convention
	}
	NextGridThisLevel = NULL; 
	NextGridNextLevel = NULL; 
	//BaryonFileName = "None"; 	

}

int grid::SetGridID(int myGridID){
	GridID = myGridID; 
	return 1; 
}

int grid::SetTask(int myTask){
	Task = myTask; 
	return 1;

}

int grid::SetGridRank(int myGridRank){
	GridRank = myGridRank; 
	return 1;
}

int grid::SetGridDimension(int myGridDimension[3]){
	for(int i = 0; i < 3; i++){
		GridDimension[i] = myGridDimension[i]; 
	}
	return 1; 
}

int grid::SetGridStartIndex(int myGridStartIndex[3]){
	for(int i = 0; i < 3; i++){
		GridStartIndex[i] = myGridStartIndex[i]; 
	}
	return 1;
}

int grid::SetGridEndIndex(int myGridEndIndex[3]){
	for(int i = 0; i < 3; i++){
		GridEndIndex[i] = myGridEndIndex[i];
	}
	return 1;
}

int grid::SetGridLeftEdge(float myGridLeftEdge[3]){
	for(int i = 0; i < 3; i++){
		GridLeftEdge[i] = myGridLeftEdge[i];
	}
	return 1;
	}

int grid::SetGridRightEdge(float myGridRightEdge[3]){
	for(int i = 0; i < 3; i++){
		GridRightEdge[i] = myGridRightEdge[i]; 
	}
	return 1;
}

int grid::SetGridLevel(int myGridLevel){
	GridLevel = myGridLevel; 
	return 1;
}

int grid::SetTime(float myTime){
	Time = myTime; 
	return 1; 
}

int grid::SetSubgridsAreStatic(int mySubgridsAreStatic){
	SubgridsAreStatic = mySubgridsAreStatic; 
	return 1; 
}

int grid::SetNumberOfBaryonFields(int myNumberOfBaryonFields){
	NumberOfBaryonFields = myNumberOfBaryonFields; 
	return 1; 
}

int grid::SetFieldType(int myFieldType[100]){
	for(int i = 0; i < 100; i++){
		FieldType[i] = myFieldType[i]; 
	}
	return 1;
}

int grid::SetBaryonFileName(char myBaryonFileName[100]){
    for(int i = 0; i < 100; i++){
	    BaryonFileName[i] = myBaryonFileName[i];
       }
	return 1;	
}

int grid::SetCourantSafetyNumber(float myCourantSafetyNumber){
	CourantSafetyNumber = myCourantSafetyNumber; 
	return 1; 
}

int grid::SetPPMFlatteningParameter(int myPPMFlatteningParameter){
	PPMFlatteningParameter = myPPMFlatteningParameter;
	return 1;
}

int grid::SetPPMDiffusionParameter(int myPPMDiffusionParameter){
	PPMDiffusionParameter = myPPMDiffusionParameter; 
	return 1; 
}

int grid::SetPPMSteepeningParameter(int myPPMSteepeningParameter){
	PPMSteepeningParameter = myPPMSteepeningParameter; 
	return 1; 
}

int grid::SetNumberOfParticles(int myNumberOfParticles){
	NumberOfParticles = myNumberOfParticles; 
	return 1; 
}

int grid::SetNumberOfActiveParticles(int myNumberOfActiveParticles){
	NumberOfActiveParticles = myNumberOfActiveParticles; 
	return 1;
}

int grid::SetGravityBoundaryType(int myGravityBoundaryType){
	GravityBoundaryType = myGravityBoundaryType; 
	return 1;
}

int grid::SetNextGridThisLevel(int *myNextGridThisLevel){
	if(*myNextGridThisLevel == 0){NextGridThisLevel = nullptr;}	
	NextGridThisLevel = myNextGridThisLevel; 
	return 1;
}

int grid::SetNextGridNextLevel(int *myNextGridNextLevel){
	if(myNextGridNextLevel == 0){NextGridNextLevel = nullptr;}
	else{NextGridNextLevel = myNextGridNextLevel;} 
	return 1;
}	

int grid::GetGridID(){return GridID;}
int grid::GetTask(){return Task;}
int grid::GetGridRank(){return GridRank;}
int* grid::GetGridDimension(){return GridDimension;}
int* grid::GetGridStartIndex(){return GridStartIndex;}
int* grid::GetGridEndIndex(){return GridEndIndex;}
float* grid::GetGridLeftEdge(){return GridLeftEdge;}
float* grid::GetGridRightEdge(){return GridRightEdge;}
int grid::GetGridLevel(){return GridLevel;}
float grid::GetTime(){return Time;}
int grid::GetSubgridsAreStatic(){return SubgridsAreStatic;}
int grid::GetNumberOfBaryonFields(){return NumberOfBaryonFields;}
int* grid::GetFieldType(){return FieldType;}
char* grid::GetBaryonFileName(){return BaryonFileName;}
float grid::GetCourantSafetyNumber(){return CourantSafetyNumber;}
int grid::GetPPMFlatteningParameter(){return PPMFlatteningParameter;}
int grid::GetPPMDiffusionParameter(){return PPMDiffusionParameter;}
int grid::GetPPMSteepeningParameter(){return PPMSteepeningParameter;}
int grid::GetNumberOfParticles(){return NumberOfParticles;}
int grid::GetNumberOfActiveParticles(){return NumberOfActiveParticles;}
int grid:: GetGravityBoundaryType(){return GravityBoundaryType;}
int* grid::GetNextGridThisLevel(){return NextGridThisLevel;}
int* grid::GetNextGridNextLevel(){return NextGridNextLevel;}

//function to print all data a grid has, useful for debugging
void grid::PrintAllData(){
	if(GridID != 0){
	std::cout << "GRID SUMMARY" << std::endl;
	std::cout << "GridID: " << GridID << std::endl; 
	std::cout << GridID << "Task: " << Task << std::endl; 
	std::cout << GridID << "GridRank: " << GridRank << std::endl; 
	std::cout << GridID << "GridDimension: " << GridDimension[0] << " " << GridDimension[1] << " " << GridDimension[2] << std::endl;
	std::cout << GridID << "GridStartIndex: " << GridStartIndex[0] << " " << GridStartIndex[1] << " " << GridStartIndex[2] << std::endl;
	std::cout << GridID << "GridEndIndex: " << GridEndIndex[0] << " " << GridEndIndex[1] << " " << GridEndIndex[2] << std::endl;
	std::cout << GridID << "GridLeftEdge: " << GridLeftEdge[0] << " " << GridLeftEdge[1] << " " << GridLeftEdge[2] << std::endl;
	std::cout << GridID << "GridRightEdge: " << GridRightEdge[0] << " " << GridRightEdge[1] << " " << GridRightEdge[2] << std::endl;
	std::cout << GridID << "GridLevel: " << GridLevel << std::endl; 
	std::cout << GridID << "Time: " << Time << std::endl; 
	std::cout << GridID << "SubgridsAreStatic: " << SubgridsAreStatic << std::endl; 
	std::cout << GridID << "NumberOfBaryonFields: " << NumberOfBaryonFields << std::endl;
	std::cout << GridID << "FieldType: "; 
        for(int i = 0; i < NumberOfBaryonFields; i++){
		std::cout << FieldType[i] << " ";
	}
	std::cout << std::endl;	
    if(BaryonFileName == NULL){
        std::cout << "BaryonFileName: NULL" << endl;
    }
    else{
	    std::cout << "BaryonFileName: " << BaryonFileName << std::endl; 
    }
	std::cout << GridID << "CourantSafetyNumber: " << CourantSafetyNumber << std::endl; 
	std::cout << GridID << "PPMFlatteningParameter: " << PPMFlatteningParameter << std::endl; 
	std::cout << GridID << "PPMDiffusionParameter: " << PPMDiffusionParameter << std::endl; 
	std::cout << GridID << "PPMSteepeningParamter: " << PPMSteepeningParameter << std::endl; 
	std::cout << GridID << "NumberOfParticles: " << NumberOfParticles << std::endl; 
	std::cout << GridID << "NumberOfActiveParticles: " << NumberOfActiveParticles << std::endl; 
	std::cout << GridID << "GravityBoundaryType: " << GravityBoundaryType << std::endl;
    if(NextGridThisLevel == NULL){
        cout << GridID << "NextGridThisLevel: NULL" << endl;
    }
    else{
	    std::cout << GridID << "NextGridThisLevel: " << *NextGridThisLevel << std::endl; 
    }
    if(NextGridNextLevel == NULL){ 
        std::cout << GridID << "NextGridNextLevel: NULL" << endl;
    }
    else{
        std::cout << GridID << "NextGridNextLevel: " << *NextGridNextLevel << std::endl; 
    } 
	}

}

int ReadListOfInts(FILE *fptr, int N, int nums[], int TestGridID){
	for(int i = 0; i < N; i++){
		if(fscanf(fptr, "%d", nums + i) != 1){
			cout << "AHHHHH SOMETHING BROKEN "<< i << N << " " << TestGridID <<  endl;
			return 0;
		}
	}
	fscanf(fptr, "\n");
	return 1;
}

//function that goes through .hierarchy file and populates grid objects with data
int parse_hierarchy(const char *filename, int lower_grid, int upper_grid, grid localgrid[], int num_grids, int nodenum){
	//declarations for IDs and data 
	int TestGridID, NextGridThisLevelID, NextGridNextLevelID, Task, GridRank, SubgridsAreStatic, NumberOfBaryonFields, PointerGridID;
	int GridStartIndex[3], GridEndIndex[3], GridDimension[3], PPMFlatteningParameter, PPMDiffusionParameter, PPMSteepeningParameter; 
	int NumberOfParticles, GravityBoundaryType; 
	float GridLeftEdge[3], GridRightEdge[3],  CourantSafetyNumber;
	float Time;
	int FieldType[100]; //change to MAX_NUMBER_OF_BARYON_FIELDS eventually
	char name[200]; 
	FILE *fptr;
	fptr = fopen(filename, "r");
	if(fptr==NULL){
		cout << "Error opening file" << endl;
	}
	else{
		cout << "Opened file: " << filename << endl;
	}
	TestGridID = -1;
    int count = 0; //count of populated grids on current node 
    while(TestGridID < num_grids){
    //read grid header
	if(fscanf(fptr,"\nGrid = %d\n", &TestGridID) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			count++;
			cout << "new grid id: " << TestGridID << " nodenum: " << nodenum << endl; 
			if(localgrid[count].SetGridID(TestGridID) != 1){
				cout << "FIRE IN TESTGRIDID SET" << endl;
			}
		}
	}
    //start parsing other lines from .hierarchy
	if(fscanf(fptr, "Task = %d\n", &Task) == 1){ 
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetTask(Task) != 1){
				cout << "FIRE IN TASK SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "GridRank = %d\n", &GridRank) == 1){ 
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridRank(GridRank) != 1){
				cout << "FIRE IN GRID RANK SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "GridDimension = %d %d %d\n", &GridDimension[0], &GridDimension[1], &GridDimension[2]) == 3){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridDimension(GridDimension) != 1){
				cout << "FIRE IN GRID DIMENSION SET" << endl;	
			}	
		}
	
	}
	//read StartIndex
	if(fscanf(fptr, "GridStartIndex = %d %d %d\n", &GridStartIndex[0], &GridStartIndex[1], &GridStartIndex[2]) == 3){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridStartIndex(GridStartIndex) != 1){
				cout << "FIRE IN GRID START INDEX SET" << endl;
			}
		}	
	} 
	//read EndIndex 
	if(fscanf(fptr, "GridEndIndex = %d %d %d\n", &GridEndIndex[0], &GridEndIndex[1], &GridEndIndex[2]) == 3){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridEndIndex(GridEndIndex) != 1){
				cout << "FIRE IN GRID END INDEX SET" << endl;
			}	
		}
	}	
	//read LeftEdge 
	if(fscanf(fptr, "GridLeftEdge = %f %f %f\n", &GridLeftEdge[0], &GridLeftEdge[1], &GridLeftEdge[2]) == 3){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridLeftEdge(GridLeftEdge) != 1){
				cout << "FIRE IN GRID LEFT EDGE SET" << endl; 
			}
		}
	}
	//read RightEdge 
	if(fscanf(fptr, "GridRightEdge = %f %f %f\n", &GridRightEdge[0], &GridRightEdge[1], &GridRightEdge[2]) == 3){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGridRightEdge(GridRightEdge) != 1){
				cout << "FIRE IN GRID RIGHT EDGE SET" << endl; 
			}
		}	
	}
	if(fscanf(fptr, "Time = %f\n", &Time) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetTime(Time) != 1){
				cout << "FIRE IN TIME SET" << endl; 
			}
		}
	}
	if(fscanf(fptr, "SubgridsAreStatic = %f\n", &SubgridsAreStatic) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetSubgridsAreStatic(SubgridsAreStatic) != 1){
				cout << "FIRE IN SUBGRIDSARESTATIC SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "NumberOfBaryonFields = %d\n", &NumberOfBaryonFields) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetNumberOfBaryonFields(NumberOfBaryonFields) != 1){
				cout << "FIRE IN NUMBEROFBARYONFIELDS SET" << endl;	
			}
		}	
	}
    //parse first part of field type line
	fscanf(fptr, "FieldType = "); 
	ReadListOfInts(fptr, NumberOfBaryonFields, FieldType, TestGridID);//read array of field IDs 
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetFieldType(FieldType) != 1){
				cout << "FIRE IN FIELDTYPE SET" << endl;
			}
		}	
	
	
	if(fscanf(fptr, "BaryonFileName = %s\n", name) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetBaryonFileName(name) != 1){
				cout << "FIRE IN BARYONFILENAME SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "CourantSafetyNumber = %f\n", &CourantSafetyNumber) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetCourantSafetyNumber(CourantSafetyNumber) != 1){
				cout << "FIRE IN COURANTSAFETYNUMBER SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "PPMFlatteningParameter = %d\n", &PPMFlatteningParameter) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetPPMFlatteningParameter(PPMFlatteningParameter) != 1){
				cout << "FIRE IN PPMFLATTENINGPARAMETER SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "PPMDiffusionParameter = %d\n", &PPMDiffusionParameter) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetPPMDiffusionParameter(PPMDiffusionParameter) != 1){
				cout << "FIRE IN PPMDIFFUSIONPARAMETER SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "PPMSteepeningParameter = %d\n", &PPMSteepeningParameter) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetPPMSteepeningParameter(PPMSteepeningParameter) != 1){
				cout << "FIRE IN PPMSTEEPENINGPARAMETER SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "NumberOfParticles = %d\n", &NumberOfParticles) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetNumberOfParticles(NumberOfParticles) != 1){
				cout << "FIRE IN NUMBEROFPARTICLES SET" << endl;
			}	
		}
	}
	int NumberOfActiveParticles; 
	if(fscanf(fptr, "NumberOfActiveParticles = %d\n", &NumberOfActiveParticles) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetNumberOfActiveParticles(NumberOfActiveParticles) != 1){
				cout << "FIRE IN NUMBEROFACTIVEPARTICLES SET" << endl;
			}
		}
	}
	fscanf(fptr, "PresentParticleTypes = \n");
	fscanf(fptr, "ParticleTypeCounts = \n");
	if(fscanf(fptr, "GravityBoundaryType = %d\n", &GravityBoundaryType) == 1){
		if(TestGridID >= lower_grid && TestGridID <= upper_grid){
			if(localgrid[count].SetGravityBoundaryType(GravityBoundaryType) != 1){
				cout << "FIRE IN GRAVITYBOUNDARYTYPE SET" << endl;
			}
		}
	}
	if(fscanf(fptr, "Pointer: Grid[%d]->NextGridThisLevel = %d\n", &PointerGridID, &NextGridThisLevelID) == 2){
		if(PointerGridID >= lower_grid && PointerGridID <= upper_grid){
			int *ptr = new int(); //this is bad practice with memory allocation 
			*ptr = NextGridThisLevelID; 
			if(localgrid[count].SetNextGridThisLevel(ptr) != 1){
				cout << "FIRE IN NEXTGRIDTHISLEVEL SET" << endl;
			}
		}	
	}
	//read pointer for NextGridNextLevel 
    //this needs a loop since at each GridID, there may be more than one NextGridNextLevel pointer
    //that needs to be set or changed 
	while(fscanf(fptr, "Pointer: Grid[%d]->NextGridNextLevel = %d\n", &PointerGridID, &NextGridNextLevelID) == 2){
		for(int i = 0; i < count; i++){
			if(localgrid[i].GetGridID() == PointerGridID){
				int *ptr = new int(); 
				*ptr = NextGridNextLevelID; 
				if(localgrid[count].SetNextGridNextLevel(ptr) != 1){
					cout << "FIRE IN NEXTGRIDNEXTLEVEL" << endl;
				}
			}

		}	

		/*if(TestGridID == localgrid[i].GetGridID()){
			int *ptr = new int(); //again bad practice, maybe fix later
			*ptr = NextGridNextLevelID;
			if(localgrid[count].SetNextGridNextLevel(ptr) != 1){
				cout << "FIRE IN NEXTGRIDNEXTLEVEL" << endl;
			}
            continue;
		}
        if(TestGridID < lower_grid || TestGridID > upper_grid){continue;}
        if(TestGridID >= lower_grid && TestGridID <= upper_grid){
            int *ptr = new int();
			*ptr = NextGridNextLevelID;
            for(int j = 0; j < num_grids; j++){
                if(TestGridID == localgrid[j].GetGridID() && localgrid[j].GetGridID() != 0){
                    cout << "localgrid rank being edited " << localgrid[j].GetGridID() << " " << i << endl; 
                    localgrid[j].SetNextGridNextLevel(ptr); 
                    continue; 
                }
            }
        }
        count2++;
        cout << "GridID: " << i << "count2: " << count2 << endl; */
	}
    }
	return(1);

}


int Distribute_Grids(grid grids[],int num_nodes, int global_rank, const char *filename){
	//determining the total number of grids 
	string line;
	ifstream file(filename);
	int node_num = global_rank / num_nodes; 
	int TestGridID = -1;
	int CurrentGridID = -1;  
    while(getline(file, line)){
            if(sscanf(line.c_str(), "\nGrid = %d\n", &TestGridID) == 1){
               if(TestGridID > CurrentGridID){
                    CurrentGridID = TestGridID;
                }
            }
    }
	int total_grids = CurrentGridID;
	int lower_grid = node_num * (total_grids / num_nodes) + 1; //lower bound of indices of grids to read on each proc 
	int upper_grid = (node_num + 1) * (total_grids / num_nodes); //upper bound of indices of grids to read on each proc
    if(total_grids % 2 == 1 and node_num == num_nodes - 1){ //catch last grid if odd number of total grids
        upper_grid++; 
    }
    cout << "lower grid num: " << lower_grid << endl; 
    cout << "upper_grid num: " << upper_grid << endl;
    //fill grid objects in range for given node 
	parse_hierarchy(filename, lower_grid, upper_grid, grids, total_grids, node_num); 
	return 1; 
}
