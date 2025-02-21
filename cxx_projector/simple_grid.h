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
#include<hdf5.h>
#include<map>
#include<math.h>
using namespace std;

class grid; 

class field
{
    private: 
        float* data; 
        const char* field_name; 
        int GridID;

    public:
        friend grid;
        field(); //constructor JACOB MAKE THE DESTRUCTOR ALREADY DONT BE CRINGE
        ~field();
        //accessor functions 
        float* GetFieldData(); 
        const char* GetFieldName(); 
        int GetFieldGridID(); 

        //setter routines 
        int SetFieldName(const char *myfieldfame); 
        int SetFieldGridID(int myfieldGridID);
        int SetFieldData(float *myfielddata, int field_size); 
        int AllocateField(int size);
        int load_data(char *filename, const char *fieldname,int GridID, int size);
};

field::field(){
    data = NULL; 
    field_name = NULL; 
    GridID = 0; 
}

field::~field(){
    data = NULL;
    delete[] data; 
}

float* field::GetFieldData(){
    if(data == NULL){
        cout << "Error: The Field you're trying to access has not been read from disk!" << endl;
    }
    return data;
}
const char* field::GetFieldName(){return field_name;}
int field::GetFieldGridID(){return GridID;}

int field::AllocateField(int size){
    data = new float[size]; 
    return 1; 
}

int field::SetFieldName(const char* myfieldname){
    field_name = myfieldname; 
    return 1; 
}

int field::SetFieldGridID(int myfieldGridID){
    GridID = myfieldGridID; 
    return 1; 
}

int field::SetFieldData(float* myfielddata, int field_size){
    for(int i = 0; i < field_size; i++){
        data[i] = myfielddata[i];
    }
    return 1; 
}

int field::load_data(char *filename, const char *fieldname, int GridID, int size){
    if(GridID == 0){return 1;}
    data = new float[size];  
    hid_t file_id, group_id, dset_id; 
    hid_t h5_status;
    char group_name[100]; 
    sprintf(group_name,"Grid%08d",GridID); 
    file_id = H5Fopen(filename, H5F_ACC_RDONLY,H5P_DEFAULT);
    group_id = H5Gopen2(file_id, group_name,H5P_DEFAULT);
    dset_id = H5Dopen(group_id, fieldname, H5P_DEFAULT); 
    h5_status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    h5_status = H5Fclose(file_id);
    return 1; 
}

    



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
		grid *NextGridThisLevel;  
		grid *NextGridNextLevel;
		float *BaryonField[20]; //pointers to BaryonFields
        int *isRefined;
        int NextGridNextLevelID; 
        int NextGridThisLevelID;
        int ProcessorNumber;
        float *CellWidth[3]; //hardcoded 3D for now 
        int field_size;
        double LengthUnits;
        double MassUnits; 
        double TimeUnits; 
        field fields[99];
        float *xyz[3]; 
        float *dxyz[3];
        int number_of_derived_fields;
        int notRefined_num; 
        int *notRefined_ind; 

	public: //member functions
		grid(); //constructor which sets everything to defaults
        ~grid();
        int make_xyz();
        int find_field_ind(int myfieldid); 

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
		int SetNextGridThisLevel(grid *myNextGridThisLevel); 
		int SetNextGridNextLevel(grid *myNextGridNextLevel);
        int SetNextGridNextLevelID(int myNextGridNextLevelID); 
        int SetNextGridThisLevelID(int myNextGridThisLevelID); 
        int SetProcessorNumber(int myProcessorNumber);
        int SetSize(int mysize);
        int SetLengthUnits(double myLengthUnits);
        int SetMassUnits(double myMassUnits); 
        int SetTimeUnits(double myTimeUnits); 
        int allocate_field_size(int size, int field_ind);
        int SetDerivedField(const char* fieldname, int (*func)(grid, float*)); 
        int allocate_isRefined(int size); 
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
		grid* GetNextGridThisLevel(); 
		grid* GetNextGridNextLevel();
	    void PrintAllData();
        int GetNextGridThisLevelID(); 
        int GetNextGridNextLevelID();
        int GetProcessorNumber();
        int GetSize();
        double GetLengthUnits();
        double GetMassUnits(); 
        double GetTimeUnits(); 
        map<string, int> field_id{{"Density", 0}, {"TotalEnergy", 1}, {"InternalEnergy", 2}, {"Pressure", 3}, {"Velocityx", 4}, {"Velocityy", 5}, {"Velocityz", 6}, {"ElectronDensity",7}, {"HIDensity",8}, {"HIIDensity",9}, {"HeIDensity", 10}, {"HeIIDensity", 11}, {"HeIIIDensity", 12}, {"HMDensity", 13}, {"H2IDensity", 14}, {"H2IIDensity", 15}, {"Metallicity",20}, {"Bx", 49}, {"By", 50}, {"Bz", 51}};
        float* GetPrimativeField(const char* field_name);  

		//friend functions
		friend int parse_hierarchy(const char *filename, grid localgrid[], int nodenum); 
		friend int ReadListOfInts(FILE *fptr, int N, int num[], int TestGridID);
		friend int distribute_grids(grid grids[],int num_nodes, int global_rank, const char *filename);
        friend int read_dataset(grid grids[], int size, int global_rank, int grid_num, const char* field_name);
        friend int isRefinedMaster(grid* current_grid);
        friend int isRefinedWorker(grid* current_grid, grid* subgrid);
        friend int FindGridID(grid grids[], int found_grids, int GridID);
        friend int SetUnits(const char *filename, grid currentgrid);
        friend int test_derived_field(grid grid, float res[]); 
        friend int create_derived_field(grid localgrids[], int global_rank, int grid_num, const char* field_name);
        //friend int SetNumNotRefined(grid mygrid); 
        //friend int SetNotRefinedInds(grid mygrid);
};


int GetNumberOfGrids(const char* filename); 

int GetNumberOfGrids(const char* filename){
	string line;
	ifstream file(filename);
	int TestGridID = -1;
	int CurrentGridID = -1;  
    while(getline(file, line)){
            if(sscanf(line.c_str(), "\nGrid = %d\n", &TestGridID) == 1){
               if(TestGridID > CurrentGridID){
                    CurrentGridID = TestGridID;
                }
            }
    }
    return CurrentGridID; 
}

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
    for(int i = 0; i < 3; i++){
        CellWidth[i] = NULL;
    }
    field_size = 0;
    number_of_derived_fields = 0;
    isRefined = NULL;
    notRefined_ind = NULL;
    notRefined_num = 0;
	//BaryonFileName = "None"; 	

}

grid::~grid(){
    for(int dim = 0; dim < GridRank; dim++){
        delete[] dxyz[dim]; 
        delete[] xyz[dim];
        CellWidth[dim] = NULL;
        delete[] CellWidth[dim]; 
    }
    isRefined = NULL;
    delete[] isRefined; 
}

int test_derived_field(grid grid, float res[]){
    float* Bx_arr = grid.GetPrimativeField("Bx"); 
    float* By_arr = grid.GetPrimativeField("By");
    int fld_size = grid.GetSize();
    for(int i = 0; i < fld_size; i++){
        //cout << "setting derived field index " << i << endl;
        //cout << "Bx[i] By[i] " << Bx_arr[i] << By_arr[i] << endl;
        //res[i] = Bx_arr[i] + By_arr[i]; 
        res[i] = 1;
    }
    return 1; 
}

float* grid::GetPrimativeField(const char* field_name){
    int fld_id = field_id[field_name];
    int field_ind = find_field_ind(fld_id);
    return fields[field_ind].GetFieldData(); 
}

int grid::find_field_ind(int myfieldid){
    int fieldid = 99; 
    for(int i = 0; i < NumberOfBaryonFields + number_of_derived_fields; i++){
        if(FieldType[i] == myfieldid){
            fieldid = i; 
        }
    }
    return fieldid; 
}

int grid::allocate_field_size(int size, int field_ind){
   fields[field_ind].data = new float[size]; 
   return 1; 
}

int grid::SetDerivedField(const char* field_name, int (*func)(grid, float*)){
   if(GridID == 0){return 1;} //Nothing to add to empty grid
   int derived_field_index = 90 + number_of_derived_fields; 
   field_id.insert({field_name, derived_field_index});
   float derived_field_data[field_size];
   func(*this, derived_field_data);
   fields[derived_field_index].AllocateField(field_size); 
   fields[derived_field_index].SetFieldData(derived_field_data, field_size);
   fields[derived_field_index].SetFieldName(field_name); 
   return 1;
}

int grid::make_xyz(){
    if(field_size == 0){
        for(int i = 0; i < GetGridRank(); i++){field_size *= GridDimension[i];}
    }
    for(int dim = 0; dim < 3; dim++){
        dxyz[dim] = new float[field_size]; 
        xyz[dim] = new float[field_size];
    }
    float GridWidth[3];
    for(int dim = 0; dim < GetGridRank(); dim++){
        GridWidth[dim] = GridRightEdge[dim] - GridLeftEdge[dim]; 
    }
    for(int k = 0; k < GridDimension[2]; k++){
        for(int j = 0; j < GridDimension[1]; j++){
            for(int i = 0; i < GridDimension[0]; i++){
                int index = (k*GridDimension[1] + j)*GridDimension[0] + i;     
                for(int dim = 0; dim < GridRank; dim++){
                    dxyz[dim][index] = GridWidth[dim]/GridDimension[dim];  
                xyz[0][index] = GridDimension[dim] * i + dxyz[index][0];  
                xyz[1][index] = GridDimension[dim] * j + dxyz[index][1];  
                xyz[2][index] = GridDimension[dim] * k + dxyz[index][2];  
                }
            }
        }
    }
    //now convert to cgs units 
    for(int dim = 0; dim < GetGridRank(); dim++){
        for(int ind = 0; ind < field_size; ind ++){
            dxyz[dim][ind] *= LengthUnits; 
            xyz[dim][ind] *= LengthUnits; 
        }
    }
    return 1; 

}

int grid::SetLengthUnits(double myLengthUnits){
    LengthUnits = myLengthUnits;
    return 1; 
}

int grid::SetMassUnits(double myMassUnits){
    MassUnits = myMassUnits; 
    return 1; 
}

int grid::SetTimeUnits(double myTimeUnits){
    TimeUnits = myTimeUnits; 
    return 1;
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

int grid::SetNextGridThisLevel(grid *myNextGridThisLevel){
	NextGridThisLevel = myNextGridThisLevel; 
	return 1;
}

int grid::SetNextGridNextLevel(grid *myNextGridNextLevel){
	NextGridNextLevel = myNextGridNextLevel;
	return 1;
}

int grid::SetNextGridThisLevelID(int myNextGridThisLevelID){
    NextGridThisLevelID = myNextGridThisLevelID;
    return 1;
}

int grid::SetNextGridNextLevelID(int myNextGridNextLevelID){
    NextGridNextLevelID = myNextGridNextLevelID;
    return 1;
}

int grid::SetProcessorNumber(int myProcessorNumber){
    ProcessorNumber = myProcessorNumber; 
    return 1;
}

int grid::SetSize(int mysize){ 
    field_size = mysize; 
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
grid* grid::GetNextGridThisLevel(){return NextGridThisLevel;}
grid* grid::GetNextGridNextLevel(){return NextGridNextLevel;}
int grid::GetNextGridThisLevelID(){return NextGridThisLevelID;}
int grid::GetNextGridNextLevelID(){return NextGridNextLevelID;}
int grid::GetProcessorNumber(){return ProcessorNumber;}
int grid::GetSize(){return field_size;}
double grid::GetLengthUnits(){return LengthUnits;}
double grid::GetMassUnits(){return MassUnits;}
double grid::GetTimeUnits(){return TimeUnits;}

//function to print all data a grid has, useful for debugging
void grid::PrintAllData(){
if(GridID > 0){
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
	    std::cout << GridID << "NextGridThisLevel: " << NextGridThisLevel->GridID << std::endl; 
    }
    if(NextGridNextLevel == NULL){ 
        std::cout << GridID << "NextGridNextLevel: NULL" << endl;
    }
    else{
        std::cout << GridID << "NextGridNextLevel: " << NextGridNextLevel->GridID << std::endl; 
    } 
    std::cout << "NextGridThisLevelID: " << NextGridThisLevelID << endl;
    std::cout << "NextGridNextLevelID: " << NextGridNextLevelID << endl;
    std::cout << "ProcessorNumber: " << ProcessorNumber << endl;
    std::cout << "size: " << field_size << endl;
    std::cout << GridID << "isRefined: "; 
    if(isRefined == NULL){cout << "NULL";}
    else{
        for(int i = 0; i < 20000; i++){
            std::cout << isRefined[i] << " ";
        }
    }
    std::cout << endl;
    std::cout << "LengthUnits " << LengthUnits << endl; 
    std::cout << "TimeUnits " << TimeUnits << endl;
    std::cout << "MassUnits " << MassUnits << endl;
    cout << "field_id " << endl;
    for(const auto& pair : field_id){
        cout << pair.first << " : " << pair.second << endl;
    }
}	

}

int grid::allocate_isRefined(int size){
    isRefined = new int[size]; 
    return 1;
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

int FindGridID(grid grids[], int found_grids, int GridID){
    int index = -1;
   for(int i = 0; i < found_grids; i++){
        if(grids[i].GetGridID() == GridID){
            index = i;
        }
    }
   return index; 
}

//function that goes through .hierarchy file and populates grid objects with data
int parse_hierarchy(const char *filename, grid localgrid[], int nodenum){
	//declarations for IDs and data 
    int num_grids = GetNumberOfGrids(filename);
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
    int count = -1; //count of populated grids on current node 
    while(TestGridID < num_grids){
    //read grid header
	if(fscanf(fptr,"\nGrid = %d\n", &TestGridID) == 1){
			count++;
			if(localgrid[count].SetGridID(TestGridID) != 1){
				cout << "FIRE IN TESTGRIDID SET" << endl;
			}
		
	}
    //start parsing other lines from .hierarchy
	if(fscanf(fptr, "Task = %d\n", &Task) == 1){ 
			if(localgrid[count].SetTask(Task) != 1){
				cout << "FIRE IN TASK SET" << endl;
			}
	}
	if(fscanf(fptr, "GridRank = %d\n", &GridRank) == 1){ 
			if(localgrid[count].SetGridRank(GridRank) != 1){
				cout << "FIRE IN GRID RANK SET" << endl;
			}
		
	}
	if(fscanf(fptr, "GridDimension = %d %d %d\n", &GridDimension[0], &GridDimension[1], &GridDimension[2]) == 3){
			if(localgrid[count].SetGridDimension(GridDimension) != 1){
				cout << "FIRE IN GRID DIMENSION SET" << endl;	
			}	
		
	
	}
	//read StartIndex
	if(fscanf(fptr, "GridStartIndex = %d %d %d\n", &GridStartIndex[0], &GridStartIndex[1], &GridStartIndex[2]) == 3){
			if(localgrid[count].SetGridStartIndex(GridStartIndex) != 1){
				cout << "FIRE IN GRID START INDEX SET" << endl;
			}
			
	} 
	//read EndIndex 
	if(fscanf(fptr, "GridEndIndex = %d %d %d\n", &GridEndIndex[0], &GridEndIndex[1], &GridEndIndex[2]) == 3){
			if(localgrid[count].SetGridEndIndex(GridEndIndex) != 1){
				cout << "FIRE IN GRID END INDEX SET" << endl;
			}	
		
	}	
	//read LeftEdge 
	if(fscanf(fptr, "GridLeftEdge = %f %f %f\n", &GridLeftEdge[0], &GridLeftEdge[1], &GridLeftEdge[2]) == 3){
			if(localgrid[count].SetGridLeftEdge(GridLeftEdge) != 1){
				cout << "FIRE IN GRID LEFT EDGE SET" << endl; 
			}
		
	}
	//read RightEdge 
	if(fscanf(fptr, "GridRightEdge = %f %f %f\n", &GridRightEdge[0], &GridRightEdge[1], &GridRightEdge[2]) == 3){
			if(localgrid[count].SetGridRightEdge(GridRightEdge) != 1){
				cout << "FIRE IN GRID RIGHT EDGE SET" << endl; 
			}
			
	}
	if(fscanf(fptr, "Time = %f\n", &Time) == 1){
			if(localgrid[count].SetTime(Time) != 1){
				cout << "FIRE IN TIME SET" << endl; 
			}
		
	}
	if(fscanf(fptr, "SubgridsAreStatic = %d\n", &SubgridsAreStatic) == 1){
			if(localgrid[count].SetSubgridsAreStatic(SubgridsAreStatic) != 1){
				cout << "FIRE IN SUBGRIDSARESTATIC SET" << endl;
			}
		
	}
	if(fscanf(fptr, "NumberOfBaryonFields = %d\n", &NumberOfBaryonFields) == 1){
			if(localgrid[count].SetNumberOfBaryonFields(NumberOfBaryonFields) != 1){
				cout << "FIRE IN NUMBEROFBARYONFIELDS SET" << endl;	
			}
			
	}
    //parse first part of field type line
	fscanf(fptr, "FieldType = "); 
	ReadListOfInts(fptr, NumberOfBaryonFields, FieldType, TestGridID);//read array of field IDs 
			if(localgrid[count].SetFieldType(FieldType) != 1){
				cout << "FIRE IN FIELDTYPE SET" << endl;
			}
			
	
	
	if(fscanf(fptr, "BaryonFileName = %s\n", name) == 1){
			if(localgrid[count].SetBaryonFileName(name) != 1){
				cout << "FIRE IN BARYONFILENAME SET" << endl;
			}
		
	}
	if(fscanf(fptr, "CourantSafetyNumber = %f\n", &CourantSafetyNumber) == 1){
			if(localgrid[count].SetCourantSafetyNumber(CourantSafetyNumber) != 1){
				cout << "FIRE IN COURANTSAFETYNUMBER SET" << endl;
			}
		
	}
	if(fscanf(fptr, "PPMFlatteningParameter = %d\n", &PPMFlatteningParameter) == 1){
			if(localgrid[count].SetPPMFlatteningParameter(PPMFlatteningParameter) != 1){
				cout << "FIRE IN PPMFLATTENINGPARAMETER SET" << endl;
			}
		
	}
	if(fscanf(fptr, "PPMDiffusionParameter = %d\n", &PPMDiffusionParameter) == 1){
			if(localgrid[count].SetPPMDiffusionParameter(PPMDiffusionParameter) != 1){
				cout << "FIRE IN PPMDIFFUSIONPARAMETER SET" << endl;
			}
		
	}
	if(fscanf(fptr, "PPMSteepeningParameter = %d\n", &PPMSteepeningParameter) == 1){
			if(localgrid[count].SetPPMSteepeningParameter(PPMSteepeningParameter) != 1){
				cout << "FIRE IN PPMSTEEPENINGPARAMETER SET" << endl;
			}
		
	}
	if(fscanf(fptr, "NumberOfParticles = %d\n", &NumberOfParticles) == 1){
			if(localgrid[count].SetNumberOfParticles(NumberOfParticles) != 1){
				cout << "FIRE IN NUMBEROFPARTICLES SET" << endl;
			}	
		
	}
	int NumberOfActiveParticles; 
	if(fscanf(fptr, "NumberOfActiveParticles = %d\n", &NumberOfActiveParticles) == 1){
			if(localgrid[count].SetNumberOfActiveParticles(NumberOfActiveParticles) != 1){
				cout << "FIRE IN NUMBEROFACTIVEPARTICLES SET" << endl;
			}
		
	}
	fscanf(fptr, "PresentParticleTypes = \n");
	fscanf(fptr, "ParticleTypeCounts = \n");
	if(fscanf(fptr, "GravityBoundaryType = %d\n", &GravityBoundaryType) == 1){
			if(localgrid[count].SetGravityBoundaryType(GravityBoundaryType) != 1){
				cout << "FIRE IN GRAVITYBOUNDARYTYPE SET" << endl;
			}
		
	}
    //read ID of NextGridThisLevel
	if(fscanf(fptr, "Pointer: Grid[%d]->NextGridThisLevel = %d\n", &PointerGridID, &NextGridThisLevelID) == 2){
            for(int i = 0; i <= count; i++){
                if(localgrid[i].GetGridID() == PointerGridID){
			if(localgrid[FindGridID(localgrid, count + 1, PointerGridID)].SetNextGridThisLevelID(NextGridThisLevelID) != 1){
				cout << "FIRE IN NEXTGRIDTHISLEVEL SET" << endl;
			}
                }
            }
    //read ID of NextGridNextLevel 
	while(fscanf(fptr, "Pointer: Grid[%d]->NextGridNextLevel = %d\n", &PointerGridID, &NextGridNextLevelID) == 2){
		for(int i = 0; i <= count; i++){
			if(localgrid[i].GetGridID() == PointerGridID){
                if(localgrid[i].SetNextGridNextLevelID(NextGridNextLevelID) != 1){
					cout << "FIRE IN NEXTGRIDNEXTLEVEL" << endl;
				}
			}

		}	

	}
}
}
    
    //Now actually set the pointers for NGTL and NGNL
    for(int i = 0; i < count; i++){
        int ind = FindGridID(localgrid, count, localgrid[i].GetNextGridThisLevelID());
        grid *ptr; 
        if(ind == -1){ptr = NULL;}
        else{ptr = localgrid + ind;}
        if(localgrid[i].SetNextGridThisLevel(ptr) != 1){
            cout << "CANT SET NGTL PTR RIGHT" << endl;
        }
    }
    for(int i = 0; i < count; i++){
        int ind = FindGridID(localgrid, count, localgrid[i].GetNextGridNextLevelID());
        grid *ptr; 
        if(ind == -1){ptr = NULL;}
        else{ptr = localgrid + ind;}
        if(localgrid[i].SetNextGridNextLevel(ptr) != 1){
            cout << "CANT SET NGNL PTR RIGHT" << endl;
        }
   }
    //set cellwidth array 
    for(int i = 0; i < count; i++){
        for(int j = 0; j < localgrid[i].GridRank; j++){
            localgrid[i].CellWidth[j] = new float[localgrid[i].GridDimension[j]];
            float delta = (localgrid[i].GridRightEdge[j] - localgrid[i].GridLeftEdge[j]) / localgrid[i].GridDimension[j];
            //put delta in array for CellWidth 
            for(int k = 0; k < localgrid[i].GridDimension[j]; k++){
                localgrid[i].CellWidth[j][k] = delta; 
            }
        }
    }

    //set units 
    for(int i = 0; i < count; i++){
        cout << "SetUnits for grid " << localgrid[i].GridID << endl;
        SetUnits("TT0000/time0000", localgrid[i]); 
        cout << "localgrid[i] LengthUnits " << localgrid[i].LengthUnits << endl;
    }

    //set field_size
    for(int i = 0; i <= count; i++){
        int size = 1; 
        for(int dim = 0; dim < localgrid[i].GridRank; dim++){
            size *= localgrid[i].GridDimension[dim]; 
        }
        localgrid[i].field_size = size; 
    }
	return(1);

}

int Distribute_Grids(grid grids[],int num_nodes, int global_rank, const char *filename){
	//determining the total number of grids 
	int node_num = global_rank / num_nodes; 
	int total_grids = GetNumberOfGrids(filename);
	int lower_grid = node_num * (total_grids / num_nodes) + 1; //lower bound of indices of grids to read on each proc 
	int upper_grid = (node_num + 1) * (total_grids / num_nodes); //upper bound of indices of grids to read on each proc
    cout << "node_num " << node_num << endl; 
    cout << "lower " << lower_grid << endl;
    cout << "upper " << upper_grid << endl;
    cout << "total grids " << total_grids << endl;
    if(total_grids % 2 == 1 and node_num == num_nodes - 1){ //catch last grid if odd number of total grids
        upper_grid++; 
    }
    //set proc num to read files from disk 
    for(int i = 0; i < total_grids; i++){
        if(grids[i].GetGridID() >= lower_grid && grids[i].GetGridID() <= upper_grid){
            grids[i].SetProcessorNumber(global_rank);
        }
    }
	return 1; 
}


int read_dataset(grid localgrids[], int num_grids, int global_rank,int grid_num, const char* field_name){
    //get length of localgrids
    //loop over localgrids
    int field_ID = 99; 
    char *filename = localgrids[grid_num].GetBaryonFileName(); 
    int *FieldType = localgrids[grid_num].GetFieldType(); 
    int *dims = localgrids[grid_num].GetGridDimension();
    int rank = localgrids[grid_num].GetGridRank();
    int GridID = localgrids[grid_num].GetGridID();
    int ProcessorNumber = localgrids[grid_num].GetProcessorNumber();
    map<string, int> field_id{{"Density", 0}, {"TotalEnergy", 1}, {"InternalEnergy", 2}, {"Pressure", 3}, {"Velocityx", 4}, {"Velocityy", 5}, {"Velocityz", 6}, {"ElectronDensity",7}, {"HIDensity",8}, {"HIIDensity",9}, {"HeIDensity", 10}, {"HeIIDensity", 11}, {"HeIIIDensity", 12}, {"HMDensity", 13}, {"H2IDensity", 14}, {"H2IIDensity", 15}, {"Metallicity",20}, {"Bx", 49}, {"By", 50}, {"Bz", 51}};
    field_ID = field_id.at(field_name); 
    if(GridID == 0){return 1;}
    if(ProcessorNumber != global_rank){return 1;} //only read on the processor number to which the grid is assigned
    //checking that all the field labels that I expect exist and were read from .hierarchy
    int field_ind = 99; 
    for(int j = 0; j < 100; j++){
        if(FieldType[j] == field_ID){
            field_ind = j; 
            break;
        }
        if(j == 99){
            cout << "CANT FIND FIELD LABEL" << endl;
            return 0;
        }
    }
    localgrids[grid_num].fields[field_ind].load_data(filename, field_name, GridID, localgrids[grid_num].field_size); 

    isRefinedMaster(localgrids);
    
    return 1;
}

int isRefinedWorker(grid* currentgrid, grid* subgrid){
    int *dims = currentgrid->GetGridDimension();
    int field_size = 1;
    //count how many zones in parent grid
    for(int i = 0; i < currentgrid->GetGridRank(); i++){field_size *= dims[i];}
    currentgrid->SetSize(field_size); //should already be set at this point
    if(currentgrid->isRefined == NULL){
        currentgrid->isRefined = new int[currentgrid->GetSize()]; //create isRefined field is it doesn't exist
        //initialize isRefined field to zero, set 1's later
        for(int i = 0; i < currentgrid->GetSize(); i++){currentgrid->isRefined[i] = 0;} 
    }
    //code mostly ripped from Grid_ZeroSolutionUnderSubgrid
    int subgrid_start[3], subgrid_end[3]; 
    for(int i = 0; i < currentgrid->GetGridRank(); i++){
        int NumberOfGhostZones = currentgrid->GridStartIndex[i]; 
        float Left = currentgrid->GridLeftEdge[i] - currentgrid->CellWidth[i][0] * NumberOfGhostZones;
        float Right = currentgrid->GridRightEdge[i] - currentgrid->CellWidth[i][currentgrid->GridDimension[i]-1] * NumberOfGhostZones;
        //cout << "Left: " << Left << " Right: " << Right << "CellWidth [i] " << currentgrid->CellWidth[i][0] << endl;
        //cout << "LeftEdge: " << currentgrid->GridLeftEdge[i] << " NumberOfGhostZones: " << NumberOfGhostZones << endl;
        //get where other grids overlap this grid and its ghost zones 
        if(subgrid->GridRightEdge[i] <= Left || subgrid->GridLeftEdge[i] >= Right){return 1;} //case of no overlap 
       subgrid_start[i] = nearbyint((subgrid->GridLeftEdge[i] - Left) / currentgrid->CellWidth[i][0]); 
       subgrid_end[i] = nearbyint((subgrid->GridRightEdge[i] - Left) / currentgrid->CellWidth[i][0]) - 1;
        //cout << "CellWidth[i][0] " << currentgrid->CellWidth[i][0] << endl;
        subgrid_start[i] = max(subgrid_start[i],0); 
        subgrid_end[i] = min(subgrid_end[i], currentgrid->GridDimension[i]-1); 
    }
    for(int k = subgrid_start[2]; k <= subgrid_end[2]; k++)
        for(int j = subgrid_start[1]; j <= subgrid_end[1]; j++)
            for(int i = subgrid_start[0]; i <= subgrid_end[0]; i++){
               int index = (k*currentgrid->GridDimension[1] + j)*currentgrid->GridDimension[0] + i;
               if(currentgrid->isRefined[index] == 0){
                   currentgrid->isRefined[index] = 1;
               }
            }
    return 1;
}

int isRefinedMaster(grid* grids){
    for(int i = 0; i < 1780; i++){
        grid *currentgrid = grids + i;
        //if NGNL doesn't exist for first grid the loop never iterates 
        //might need to call this differently in read_grids
        grid *temp = currentgrid->NextGridNextLevel; //pointer to NextGridNextLevel
        while(temp != NULL){
            isRefinedWorker(currentgrid, temp);//Set subgrid field 
            currentgrid = temp;
            temp = currentgrid->NextGridNextLevel;//move down one in the hierarchy               
        }
    }
    return 1;
}


int SetUnits(const char* filename, grid currentgrid){
    FILE *fptr; 
    double length_units, mass_units, time_units;
    fptr = fopen(filename, "r"); 
    if(fptr==NULL){
        cout << "Error opening Units file" << endl;
        return 0; 
    }
    else{
        char line[256]; 
        while(fgets(line, sizeof(line),fptr)){
            if(sscanf(line, "LengthUnits = %lf", &length_units) == 1){
                currentgrid.SetLengthUnits(length_units);
                cout << "CurrentGridLengthUnits: " << currentgrid.LengthUnits << endl;
            }
            if(sscanf(line, "MassUnits = %lf", &mass_units) == 1){
                currentgrid.SetMassUnits(mass_units);
            }
            if(sscanf(line, "TimeUnits = %lf", &time_units) == 1){
                currentgrid.SetTimeUnits(time_units);
            }
        }
    }
    fclose(fptr);
    return 1; 
}
/*
int SetNumNotRefined(grid *mygrid){
    for(int i = 0; i < mygrid->GetSize(); i++){
        if(mygrid->isRefined[i] == 0){mygrid->notRefined_num++;}
    }
    return 1; 
}

int SetNotRefinedInds(grid *mygrid){
    mygrid->notRefined_ind = new int[mygrid->notRefined_num];
    int ind = 0;
    for(int i = 0; i < mygrid->GetSize(); i++){
        if(mygrid->isRefined[i] == 0){
            notRefined_ind[ind] = i; 
            ind++;
        }
    }
    return 1;
}*/
