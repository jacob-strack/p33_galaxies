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
#include<cstdarg>
#include<vector>
#include <png.h>

using namespace std;
int sign(float num); 
int nint(float num);

class grid; 
int toggle = 0;  //0 does isRefined with ActiveSize, 1 does full field_size with flagging out ghost zones
                 //in theory, this shouldn't change how anything runs 
class flat_array
{
    private: 
        double* data; 
        const char* field_name;
        int fsize; 

    public:
        friend grid;
        flat_array();
        flat_array(int size);
        ~flat_array();
        //accessor functions 
        double* GetFieldData(); 
        const char* GetFieldName();
        int GetSize();
        double GetDataAtInd(int ind);
        double GetTotal(); 

        flat_array operator=(const flat_array& obj){
            if(this == &obj){return *this;}
            fsize = obj.fsize; 
            field_name = obj.field_name;
            cout << "Tried Setting new field_name" << endl;
            data = new double[fsize];
            for(int i = 0; i < fsize; i ++){
                data[i] = obj.data[i]; 
            }
            return *this;
        }

        flat_array operator*(const flat_array& obj) const {
            flat_array ans(this->fsize);
            if(this->fsize != obj.fsize){return 0;}
            for(int i = 0; i < this->fsize; i++){
               ans.data[i] = this->data[i] * obj.data[i];  
            } 
            ans.fsize = this->fsize; 
            return ans; 
        }

        flat_array operator/(const flat_array& obj) const {
            flat_array ans(this->fsize);
            if(this->fsize != obj.fsize){return 0;}
            for(int i = 0; i < this->fsize; i++){
               ans.data[i] = this->data[i] / obj.data[i];  
            } 
            ans.fsize = this->fsize; 
            return ans; 
        }
        
        flat_array operator+(const flat_array& obj) const {
            flat_array ans(this->fsize);
            if(this->fsize != obj.fsize){return 0;}
            for(int i = 0; i < this->fsize; i++){
               ans.data[i] = this->data[i] + obj.data[i];  
            } 
            ans.fsize = this->fsize; 
            return ans; 
        }
        
        flat_array operator-(const flat_array& obj) const {
            flat_array ans(this->fsize);
            if(this->fsize != obj.fsize){return 0;}
            for(int i = 0; i < this->fsize; i++){
               ans.data[i] = this->data[i] - obj.data[i];  
            } 
            ans.fsize = this->fsize; 
            return ans; 
        }
        //setter routines 
        int SetFieldName(const char *myfieldfame); 
        int SetFieldData(double *myfielddata, int field_size); 
        int AllocateField(int mysize);
        int load_data(char *filename, const char *fieldname,int GridID, int size);
        int SetDataAtInd(float data_pt, int ind);
        int SetPrimative(grid *grids, int num_grids, const char *fieldname);
        int Makexyz(grid *grids, int num_grids, int dim);
        int Makedxyz(grid *grids, int num_grids, int dim);
        int MakeCellVolume(flat_array *dx, flat_array *dy, flat_array *dz);
        int SetDerivedFlatArray(double (*func)(int,double*,va_list), double* in1, ...);
        int WriteData(const char *filename); 
        
};

flat_array::flat_array(){
    data = NULL; 
    field_name = NULL;
    fsize = 0;
}
const char* flat_array::GetFieldName(){return field_name;}
flat_array::flat_array(int mysize){
    data = new double[mysize]; 
    field_name = NULL; 
    fsize = mysize; 
}


flat_array::~flat_array(){
    data = NULL;
    delete[] data; 
}
int flat_array::GetSize(){return fsize;}

int flat_array::WriteData(const char *filename){
    hid_t file_id, dset_id; 
    hid_t float_type_id, file_type_id, file_dsp_id; 
    const int rank = 1;
    const long unsigned int size = (long unsigned int)fsize;
    hsize_t OutDims[rank] = {size};
    herr_t h5_status; 
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    file_dsp_id = H5Screate_simple(1, OutDims, NULL); 
    dset_id = H5Dcreate2(file_id, "flat_array", H5T_NATIVE_FLOAT, file_dsp_id, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
    //had float_type_id next line replaced with H5T_NATIVE_FLOAT
    h5_status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    h5_status = H5Sclose(file_dsp_id); 
    h5_status = H5Dclose(dset_id);
    return 1;
}

double flat_array::GetTotal(){
    double total = 0.0;
    cout << "GetTotal fsize " << fsize << endl;
    for(int i = 0; i < fsize; i++){
        if(data[i] != data[i]){
            cout << "Somethin's nan " << i << endl;
        }
        total += data[i]; 
    }
    return total;
}
double flat_array::GetDataAtInd(int ind){
    if(ind > fsize){cout << "Bad Index for FA!" << endl;}
    return data[ind];}

int flat_array::SetDataAtInd(float data_pt, int ind){ 
   if(ind >= fsize){
        cout << "Index out of range! Nothing Set . . ." << endl;
        return 0; 
    }
   data[ind] = data_pt; 
   return 1; 
}

double* flat_array::GetFieldData(){
    if(data == NULL){
        cout << "Error: The Field you're trying to access has not been read from disk!" << endl;
    }
    return data;
}

int flat_array::AllocateField(int size){
    data = new double[size]; 
    return 1; 
}

int flat_array::SetFieldName(const char* myfieldname){
    field_name = myfieldname; 
    return 1; 
}

int flat_array::SetFieldData(double* myfielddata, int field_size){
    for(int i = 0; i < field_size; i++){
        data[i] = myfielddata[i];
    }
    return 1; 
}

int flat_array::load_data(char *filename, const char *fieldname, int GridID, int size){
    if(GridID == 0){return 1;}
    data = new double[size];  
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

int plot_array(flat_array x_arr, flat_array y_arr, flat_array fill_arr, int num_bins_x, int num_bins_y, float x_lower, float x_upper, float y_lower, float y_upper, int log);
    

class grid
{
	private: //member data 
		int GridID; 
		int Task;
		int GridRank; 
		int GridDimension[3]; 
		int GridStartIndex[3];
		int GridEndIndex[3]; 
		double GridLeftEdge[3]; 
		double GridRightEdge[3]; 
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
        double *xyz[3]; 
        double *dxyz[3];
        int number_of_derived_fields;
        int notRefined_num; 
        int *notRefined_ind;
        int ActiveSize;

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
		int SetGridLeftEdge(double myGridLeftEdge[3]); 
		int SetGridRightEdge(double myGridRightEdge[3]);
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
        int SetActiveSize(int myActiveSize); 
		//Accessor Routines
		int GetGridID();
        double** Get_dxyz();
	    int GetTask();	
		int GetGridRank(); 
		int* GetGridDimension();	
		int* GetGridStartIndex(); 
		int* GetGridEndIndex(); 
		double* GetGridLeftEdge(); 
		double* GetGridRightEdge(); 
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
        void showmepointerids();
        int GetNextGridThisLevelID(); 
        int GetNextGridNextLevelID();
        int GetProcessorNumber();
        int GetSize();
        double GetLengthUnits();
        double GetMassUnits(); 
        double GetTimeUnits();
        int GetActiveSize(); 
        map<string, int> field_id{{"Density", 0}, {"TotalEnergy", 1}, {"InternalEnergy", 2}, {"Pressure", 3}, {"Velocityx", 4}, {"Velocityy", 5}, {"Velocityz", 6}, {"ElectronDensity",7}, {"HIDensity",8}, {"HIIDensity",9}, {"HeIDensity", 10}, {"HeIIDensity", 11}, {"HeIIIDensity", 12}, {"HMDensity", 13}, {"H2IDensity", 14}, {"H2IIDensity", 15}, {"Metallicity",20}, {"Bx", 49}, {"By", 50}, {"Bz", 51}};
        //float* GetPrimativeField(const char* field_name); 
        int GetNumNotRefined();
        int* GetNotRefinedInd();
        double** Get_xyz();

		//friend functions
		friend int parse_hierarchy(const char *filename, grid localgrid[], int nodenum); 
		friend int ReadListOfInts(FILE *fptr, int N, int num[], int TestGridID);
		friend int distribute_grids(grid grids[],int num_nodes, int global_rank, const char *filename);
        friend int read_dataset(grid grids[], int size, int global_rank, int grid_num, const char* field_name);
        friend int isRefinedMaster(grid* current_grid);
        friend int isRefinedWorker(grid* current_grid, grid* subgrid);
        friend int FindGridID(grid grids[], int found_grids, int GridID);
        friend int SetUnits(const char *filename, grid *currentgrid);
        //friend int test_derived_field(grid grid, float res[]); 
        friend int create_derived_field(grid localgrids[], int global_rank, int grid_num, const char* field_name);
        friend int SetNumNotRefined(grid *mygrid); 
        friend int SetNotRefinedInds(grid *mygrid);
        friend int build_isRefined(grid *localgrid);
        friend int set_not_refined(grid *localgrid);
    
};

int GetNumNotRefinedGrids(grid *grids, int num_grids); 

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
    LengthUnits = 0;
    ActiveSize = 0;
    GridLevel = 0; 
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

/*int test_derived_field(grid grid, float res[]){
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
}*/

/*float* grid::GetPrimativeField(const char* field_name){
    int fld_id = field_id[field_name];
    int field_ind = find_field_ind(fld_id);
    return fields[field_ind].GetFieldData(); 
}*/

int grid::find_field_ind(int myfieldid){
    int fieldid = 99; 
    for(int i = 0; i < NumberOfBaryonFields + number_of_derived_fields; i++){
        if(FieldType[i] == myfieldid){
            fieldid = i; 
        }
    }
    return fieldid; 
}


int grid::make_xyz(){
    int ActiveDims[3]; 
    int NumberOfGhostZones = GridStartIndex[0];
    for(int i = 0; i < GetGridRank(); i++){
        ActiveDims[i] = GridEndIndex[i] - GridStartIndex[i] + 1;
    }
    if(field_size == 0){
        field_size = 1;
        for(int i = 0; i < GetGridRank(); i++){field_size *= ActiveDims[i];}
    }
    for(int dim = 0; dim < 3; dim++){
        dxyz[dim] = new double[field_size]; 
        xyz[dim] = new double[field_size];
    }
    double GridWidth[3];
    for(int dim = 0; dim < GetGridRank(); dim++){
        GridWidth[dim] = GridRightEdge[dim] - GridLeftEdge[dim]; 
    }
    for(int k = 0; k < ActiveDims[2]; k++){
        for(int j = 0; j < ActiveDims[1]; j++){
            for(int i = 0; i < ActiveDims[0]; i++){
                int index = (k*ActiveDims[1] + j)*ActiveDims[0] + i;    
                for(int dim = 0; dim < GridRank; dim++){
                    dxyz[dim][index] = GridWidth[dim]/ActiveDims[dim]; 
                }
                //There should be left edge + width here
                xyz[0][index] = GridLeftEdge[0] +  (i + 0.5) * dxyz[0][index];  
                xyz[1][index] = GridLeftEdge[1] +  (j + 0.5) * dxyz[1][index];  
                xyz[2][index] = GridLeftEdge[2] +  (k + 0.5) * dxyz[2][index];  
                
            }
        }
    }
    //now convert to cgs units 
    for(int dim = 0; dim < GetGridRank(); dim++){
        for(int ind = 0; ind < field_size; ind ++){
            if(LengthUnits == 0){cout << "LENGTH UNITS ZERO!" << " " << GridID << endl;}
            if(LengthUnits != 0){
            //dxyz[dim][ind] *= LengthUnits; 
            //xyz[dim][ind] *= LengthUnits; 
            }
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

int grid::SetActiveSize(int myActiveSize){
    ActiveSize = myActiveSize; 
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

int grid::SetGridLeftEdge(double myGridLeftEdge[3]){
	for(int i = 0; i < 3; i++){
		GridLeftEdge[i] = myGridLeftEdge[i];
	}
	return 1;
	}

int grid::SetGridRightEdge(double myGridRightEdge[3]){
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
double* grid::GetGridLeftEdge(){return GridLeftEdge;}
double* grid::GetGridRightEdge(){return GridRightEdge;}
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
int grid::GetNumNotRefined(){return notRefined_num;}
int* grid::GetNotRefinedInd(){return notRefined_ind;}
double** grid::Get_xyz(){return xyz;}
double** grid::Get_dxyz(){return dxyz;}
int grid::GetActiveSize(){return ActiveSize;}

int GetTotalNotRefined(grid *grids, int num_grids);
void grid::showmepointerids(){
    cout << "GridID: " << GridID << endl; 
    std::cout << "NextGridThisLevelID: " << NextGridThisLevelID << endl;
    std::cout << "NextGridNextLevelID: " << NextGridNextLevelID << endl;
}
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
        for(int i = 0; i < 2; i++){
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
    cout << "numNotRefined: " << notRefined_num << endl;
    cout << "notRefined_inds: ";
    cout << notRefined_ind << endl;
    cout << endl;
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
            if(GridID==1780){cout << "FindGridID " << GridID << " " << index << endl;}
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
	double GridLeftEdge[3], GridRightEdge[3],  CourantSafetyNumber;
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
	if(fscanf(fptr, "GridLeftEdge = %lf %lf %lf\n", &GridLeftEdge[0], &GridLeftEdge[1], &GridLeftEdge[2]) == 3){
			if(localgrid[count].SetGridLeftEdge(GridLeftEdge) != 1){
				cout << "FIRE IN GRID LEFT EDGE SET" << endl; 
			}
		
	}
	//read RightEdge 
	if(fscanf(fptr, "GridRightEdge = %lf %lf %lf\n", &GridRightEdge[0], &GridRightEdge[1], &GridRightEdge[2]) == 3){
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
	if(fscanf(fptr, "CourantSafetyNumber = %lf\n", &CourantSafetyNumber) == 1){
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
    for(int i = 0; i <= count; i++){
        int ind = FindGridID(localgrid, count+1, localgrid[i].GetNextGridThisLevelID());
        grid *ptr; 
        if(ind == -1){ptr = NULL;}
        else{ptr = localgrid + ind;}
        if(localgrid[i].SetNextGridThisLevel(ptr) != 1){
            cout << "CANT SET NGTL PTR RIGHT" << endl;
        }
    }
    for(int i = 0; i <= count; i++){
        int ind = FindGridID(localgrid, count+1, localgrid[i].GetNextGridNextLevelID());
        grid *ptr; 
        if(ind == -1){ptr = NULL;}
        else{ptr = localgrid + ind;}
        if(localgrid[i].SetNextGridNextLevel(ptr) != 1){
            cout << "CANT SET NGNL PTR RIGHT" << endl;
        }
   }
    //set cellwidth array 
    for(int i = 0; i <= count; i++){
        int ActiveDims[3];
        int ActiveSize = 1;
        for(int dim = 0; dim < localgrid[i].GridRank; dim++){
            ActiveDims[dim] = localgrid[i].GridDimension[dim] - 2*localgrid[i].GridStartIndex[dim]; 
            ActiveSize *= ActiveDims[dim]; 
        }
        localgrid[i].SetActiveSize(ActiveSize);
        for(int j = 0; j < localgrid[i].GridRank; j++){
            localgrid[i].CellWidth[j] = new float[localgrid[i].GridDimension[j]];
            float delta = (localgrid[i].GridRightEdge[j] - localgrid[i].GridLeftEdge[j]) / (ActiveDims[j]);
            //put delta in array for CellWidth 
            for(int k = 0; k < localgrid[i].GridDimension[j]; k++){
                localgrid[i].CellWidth[j][k] = delta; 
            }
        }
    }

    //set units 
    for(int i = 0; i <= count; i++){
        SetUnits("TT0000/time0000", localgrid + i);
    }

    //set field_size
    for(int i = 0; i <= count; i++){
        int size = 1; 
        for(int dim = 0; dim < localgrid[i].GridRank; dim++){
            size *= localgrid[i].GridDimension[dim]; 
        }
        localgrid[i].field_size = size;
        if(toggle == 1){
            localgrid[i].isRefined = new int[size]; 
            for(int l = 0; l < size; l++)
                localgrid[i].isRefined[l] = 0;
    }
        if(toggle == 0){ 
            localgrid[i].isRefined = new int[localgrid[i].ActiveSize];
            for(int l = 0; l < localgrid[i].ActiveSize; l++)
                localgrid[i].isRefined[l] = 0;

        }
    }
    //flag ghost zones in isRefined for case where isRefined is the same size as field_size 
    if(toggle == 1){
    for(int l = 0; l <= count; l++)
    for(int k = 0; k < localgrid[l].GridDimension[2]; k++)
        for(int j = 0; j <localgrid[l].GridDimension[1]; j++)
            for(int i = 0; i < localgrid[l].GridDimension[0]; i++){
                int index = (k*localgrid[l].GridDimension[1] + j)*localgrid[l].GridDimension[0] + i; 
                if(k < localgrid[l].GridStartIndex[2] || k > localgrid[l].GridEndIndex[2]){
                    localgrid[l].isRefined[index] = 2;
                }
                if(j < localgrid[l].GridStartIndex[1] || j > localgrid[l].GridEndIndex[1])
                    localgrid[l].isRefined[index] = 2; 
                if(i < localgrid[l].GridStartIndex[0] || i > localgrid[l].GridEndIndex[0])
                    localgrid[l].isRefined[index] = 2; 
    }
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

int isRefinedWorker(grid* currentgrid, grid* subgrid){
    int *dims = currentgrid->GetGridDimension();
    int ActiveDims[currentgrid->GetGridRank()];
    int ActiveSize = 1;
    for(int dim = 0; dim < currentgrid->GetGridRank(); dim++){
        ActiveDims[dim] = currentgrid->GridDimension[dim] - 2*currentgrid->GridStartIndex[dim];
        ActiveSize *= ActiveDims[dim];
    }
    if(currentgrid->isRefined == NULL){
        currentgrid->isRefined = new int[ActiveSize]; //create isRefined field is it doesn't exist
        //initialize isRefined field to zero, set 1's later
        cout << "making new isRefined . . ." << endl;
        for(int i = 0; i < ActiveSize; i++){currentgrid->isRefined[i] = 0;}
    }
    int field_size = 1;
    //count how many zones in parent grid
    for(int i = 0; i < currentgrid->GetGridRank(); i++){field_size *= (dims[i]);}
    currentgrid->SetSize(field_size); //should already be set at this point
    if(toggle == 0){
    //code mostly ripped from Grid_ZeroSolutionUnderSubgrid
    float subgrid_start[3], subgrid_end[3]; 
    for(int i = 0; i < currentgrid->GetGridRank(); i++){
        double Left = currentgrid->GridLeftEdge[i];
        double Right = currentgrid->GridRightEdge[i];
        //get where other grids overlap this grid and its ghost zones 
        if(subgrid->GridRightEdge[i] <= Left || subgrid->GridLeftEdge[i] >= Right){
            cout << "No overlap. Exiting." << endl;
            return 1;} //case of no overlap 
       subgrid_start[i] = (subgrid->GridLeftEdge[i] - Left) / currentgrid->CellWidth[i][0];
       subgrid_end[i] = ((subgrid->GridRightEdge[i] - Left) / currentgrid->CellWidth[i][0]) - 1;
       subgrid_start[i] = max(nint(subgrid_start[i]),0); 
       subgrid_end[i] = min(nint(subgrid_end[i]), ActiveDims[i] - 1);
    }
    for(int k = subgrid_start[2]; k <= subgrid_end[2]; k++)
        for(int j = subgrid_start[1]; j <= subgrid_end[1]; j++){
            for(int i = subgrid_start[0]; i <= subgrid_end[0]; i++){
                int index = (k*ActiveDims[1] + j)*ActiveDims[0] + i;
                currentgrid->isRefined[index] = 1;
            }
        }
    
    cout.flush();
    }
    if(toggle == 1){
    cout << currentgrid->GridID << " " << subgrid->GridID << endl;
    int subgrid_start[3], subgrid_end[3];
    for(int dim = 0; dim < currentgrid->GetGridRank(); dim++){
        subgrid_start[dim] = 0; 
        subgrid_end[dim] = 0; 
    }
    for(int i = 0; i < currentgrid->GetGridRank(); i++){
        //shorthand for left and right edges of parent grid
        double Left = currentgrid->GridLeftEdge[i]; 
        double Right = currentgrid->GridRightEdge[i];
        //check if the subgrid isn't under the current parent grid. I don't this this ever happens
        if(subgrid->GridRightEdge[i] <= Left || subgrid->GridLeftEdge[i] >= Right){
            cout << "No overlap. Exiting." << endl;
            return 1;} 
        //calculate the index of the parent grid where the subgrid starts and ends
        subgrid_start[i] = nint((subgrid->GridLeftEdge[i] - Left) / currentgrid->CellWidth[i][0]) + currentgrid->GridStartIndex[i];
        subgrid_end[i] = nint((subgrid->GridRightEdge[i] - Left) / currentgrid->CellWidth[i][0]) - 1 + currentgrid->GridStartIndex[i];
        //if subgrid starts before or ends after parent, set start/end to parent boundary w/o ghost zones
        subgrid_start[i] = max((subgrid_start[i]), currentgrid->GridStartIndex[i]);
        subgrid_end[i] = min((subgrid_end[i]), currentgrid->GridEndIndex[i]);
    } 
    //cout << currentgrid->GridID << " " << subgrid_end[0] - subgrid_start[0] << endl;
    //set all points above subgrid to 1
    for(int k = subgrid_start[2]; k <= subgrid_end[2]; k++)
        for(int j = subgrid_start[1]; j <= subgrid_end[1]; j++){
            for(int i = subgrid_start[0]; i <= subgrid_end[0]; i++){
                int index = (k*currentgrid->GridDimension[1] + j)*currentgrid->GridDimension[0] + i; 
                currentgrid->isRefined[index] = 1;
            }
        }
    }
    cout.flush();  
    return 1;
}

int sign(float num){
    if(num < 0){return -1;}
    if(num > 0){return 1;}
    else{return 0;}
}

int nint(float num){
    return (int)(num + 0.5*sign(num)); 
}

int isRefinedMaster(grid* currentgrid){
        //if NGNL doesn't exist for first grid the loop never iterates 
        //might need to call this differently in read_grids
        grid *temp = currentgrid->NextGridNextLevel; //pointer to NextGridNextLevel
        while(temp != NULL){
            isRefinedWorker(currentgrid, temp);//Set subgrid field 
            temp = temp->NextGridThisLevel;//move over one in the hierarchy               
        }
    return 1;
}


int SetUnits(const char* filename, grid *currentgrid){
    FILE *fptr; 
    fptr = fopen(filename, "r"); 
    if(fptr==NULL){
        cout << "Error opening Units file" << endl;
        return 0; 
    }
    else{
        char line[256]; 
        while(fgets(line, sizeof(line),fptr)){
            double length_units = 0;
            double mass_units = 0; 
            double time_units = 0;
            if(sscanf(line, "LengthUnits = %lf", &length_units) == 1){
                currentgrid->SetLengthUnits(length_units);
            }
            if(sscanf(line, "MassUnits = %lf", &mass_units) == 1){
                currentgrid->SetMassUnits(mass_units);
            }
            if(sscanf(line, "TimeUnits = %lf", &time_units) == 1){
                currentgrid->SetTimeUnits(time_units);
            }
        }
    }
    fclose(fptr);
    return 1; 
}

int SetNumNotRefined(grid *mygrid){
    //initialize various counts
    int count = 0;
    int ghost_count = 0;
    int refined_count = 0; 
    int ActiveSize  = 1;
    //if not a populated grid, do nothing
    if(mygrid->GridID == 0){
        mygrid->notRefined_num = 0; 
        return 1; 
    }
    for(int i = 0; i < mygrid->GetGridRank(); i++){
        ActiveSize *= (mygrid->GridDimension[i] - 2*mygrid->GridStartIndex[i]);
    }
    //count zones
    if(toggle == 0){
    for(int i = 0; i < ActiveSize; i++){
        if(mygrid->isRefined[i] == 0){count++;}
        if(mygrid->isRefined[i] == 1){refined_count++;}
    }
    }
    if(toggle == 1){
    for(int i = 0; i < mygrid->field_size; i++){
        if(mygrid->isRefined[i] == 0){count++;}
        if(mygrid->isRefined[i] == 1){refined_count++;}
        if(mygrid->isRefined[i] == 2){ghost_count++;}
    }
    }
    
    if(mygrid->ActiveSize != count + refined_count){cout << "ZONE FLAGGING PROBLEM" << endl;}
    mygrid->notRefined_num = count;
    return 1; 
}

int SetNotRefinedInds(grid *mygrid){
    if(mygrid->notRefined_num == 0){
        return 1; //nothing to do
    }
    if(mygrid->notRefined_ind != NULL){
        cout << "changing notRefined_ind array that already exists!" << endl;
        mygrid->notRefined_ind = NULL; 
        delete[] mygrid->notRefined_ind;
        mygrid->notRefined_ind = new int[mygrid->notRefined_num];
    }
    else{mygrid->notRefined_ind = new int[mygrid->notRefined_num];}
    int ind = 0;
    //loop through entire isRefined and store indices of unrefined cells of current grid
    //0 toggle is for using only ActiveSize array
    //1 toggle uses field_size array with ghost zones flagged to 2
    if(toggle == 0){
    for(int i = 0; i < mygrid->GetActiveSize(); i++){
        if(mygrid->isRefined[i] == 0){
            mygrid->notRefined_ind[ind++] = i; 
        }
    }

    }
    if(toggle == 1){
    //this would need to include some sort of coversion to index in active size. toggle 0 does it better
    for(int i = 0; i < mygrid->field_size; i++){
        if(mygrid->isRefined[i] == 0){
            mygrid->notRefined_ind[ind] = i; 
            ind++;
        }
    }
    }
    return 1;
}


int build_isRefined(grid *localgrid){
    isRefinedMaster(localgrid);
    return 1;
}

int set_not_refined(grid *localgrid){
    SetNumNotRefined(localgrid); 
    SetNotRefinedInds(localgrid); 
    return 1;
}

int flat_array::SetPrimative(grid *grids, int num_grids, const char *fieldname){
    int ind = 0;
    int thing = GetTotalNotRefined(grids, num_grids);
    cout << "fsize: " << thing << endl;
    data = new double[thing];
    fsize = thing;
    for(int i = 0; i < num_grids; i++){
        grid *currentgrid = grids + i;
        if(currentgrid->GetGridID() == 0){continue;}
        const char* filename = currentgrid->GetBaryonFileName();
        if(currentgrid->field_id.count(fieldname) == 0){
            cout << "Field Doesn't Exist!" << endl; 
            return 0;
        }
        float h5_data[currentgrid->GetActiveSize()];
        hid_t file_id, group_id, dset_id; 
        hid_t h5_status;
        char group_name[100]; 
        sprintf(group_name,"Grid%08d",currentgrid->GetGridID()); 
        file_id = H5Fopen(filename, H5F_ACC_RDONLY,H5P_DEFAULT);
        group_id = H5Gopen2(file_id, group_name,H5P_DEFAULT);
        dset_id = H5Dopen(group_id, fieldname, H5P_DEFAULT); 
        h5_status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, h5_data);
        h5_status = H5Fclose(file_id);
        int *inds = currentgrid->GetNotRefinedInd(); 
        for(int j = 0; j < currentgrid->GetNumNotRefined(); j++){
            if(inds[j] > currentgrid->GetActiveSize()){cout << "Reading bad things" << endl;}
            data[ind++] = float(h5_data[inds[j]]);
        }
    }
    return 1; 
}

int GetTotalNotRefined(grid *grids, int num_grids){
    int tot_notRefined = 0;
    int tot_zones = 0; 
    for(int i = 0; i < num_grids; i++){
        grid *currentgrid = grids + i; 
        tot_notRefined += currentgrid->GetNumNotRefined();
        tot_zones += (currentgrid->GetActiveSize() - currentgrid->GetNumNotRefined()); 
    }
    cout << "tot not refined tot zones " << tot_notRefined << " " << tot_zones << endl;
    return tot_notRefined; 
    

}

int flat_array::Makexyz(grid *grids, int num_grids, int dim){
    int ind = 0; 
    for(int i = 0; i < num_grids; i++){
        grid *currentgrid = grids + i; 
        if(currentgrid->GetGridID() == 0){continue;}
        currentgrid->make_xyz();
        int *NotRefinedInds = currentgrid->GetNotRefinedInd();
        double **xyz_arr = currentgrid->Get_xyz(); 
        for(int j = 0; j < currentgrid->GetNumNotRefined(); j++){
            data[ind++] = xyz_arr[dim][NotRefinedInds[j]];
        }
    }
    return 1;
}
int flat_array::Makedxyz(grid *grids, int num_grids, int dim){
    int ind = 0;
    cout << "Makedxyz size " << fsize << endl;
    for(int i = 0; i < num_grids; i++){
        grid *currentgrid = grids + i; 
        if(currentgrid->GetGridID() == 0){continue;}
        currentgrid->make_xyz();
        int *NotRefinedInds = currentgrid->GetNotRefinedInd();
        double **dxyz_arr = currentgrid->Get_dxyz();
        for(int j = 0; j < currentgrid->GetNumNotRefined(); j++){
            data[ind++] = dxyz_arr[dim][NotRefinedInds[j]];
        }
    }
    cout << "Makedxyz size2 " << fsize << endl;
    return 1;
}

int flat_array::MakeCellVolume(flat_array *dx, flat_array *dy, flat_array *dz){
    int rank = 3; 
    if(data == NULL){
        data = new double[dx->fsize]; 
    }
    cout << "dx fsize CV " << dx->fsize << endl;
    for(int j = 0; j < dx->fsize; j++){
            data[j] = dx->data[j]*dy->data[j]*dz->data[j]; 
    }
    cout << "Cell Volume Test " << data[100] << endl;
    return 1;
}
    

int flat_array::SetDerivedFlatArray(double (*func)(int, double*, va_list), double* in1, ...){
    for(int i = 0; i < fsize; i++){
        va_list args;
        va_start(args, in1);
        double ans; 
        ans = func(i, in1, args);
        data[i] = ans; 
        va_end(args);
    }
    cout << "Done setting data" << endl;
    return 1;
}

int plot_array(flat_array x_arr, flat_array y_arr, flat_array fill_arr,  int num_bins_x, int num_bins_y, float x_lower, float x_upper, float y_lower, float y_upper, int log){
    cout << "starting plot_array" << endl;
    cout << "sizes " << x_arr.GetSize() << " " << y_arr.GetSize() << " " << fill_arr.GetSize() << endl;
    float **histo = new float*[num_bins_y]; 
    for(int i = 0; i < num_bins_y; i++)
        histo[i] = new float[num_bins_x];
    int **counts = new int*[num_bins_y];
    for(int i = 0; i < num_bins_y; i++)
        counts[i] = new int[num_bins_x];
    cout << "After declaring arrays" << endl;
    float min_x;
    float min_y;
    float max_x; 
    float max_y;
    float min_data = 0; 
    float max_data = 1.0;
    double *x_data = x_arr.GetFieldData(); 
    double *y_data = y_arr.GetFieldData(); 
    double *fill_data = fill_arr.GetFieldData();
    cout << "zeroing histograms" << endl;
    for(int i = 0; i < num_bins_y; i++){
        for(int j = 0; j < num_bins_x; j++){
            histo[i][j] = 0;
            counts[i][j] = 0;
        }   
    }
    min_x = x_lower; 
    max_x = x_upper; 
    min_y = y_lower; 
    max_y = y_upper;
    cout << min_x << " " << max_x << " " << min_y << " " << max_y << endl;
    for(int i = 0; i < x_arr.GetSize(); i++){
        x_data[i] -= min_x;
        y_data[i] -= min_y;
        x_data[i] /= (max_x - min_x); 
        y_data[i] /= (max_y - min_y);
    }
    int testx = 0; 
    int testy = 0;
    for(int i = 0; i < fill_arr.GetSize(); i++){
        double x = x_data[i];
        double y = y_data[i];  
        double data = fill_data[i];
        int x_bin =  static_cast<int>(x * (num_bins_x - 1)); 
        int y_bin =  static_cast<int>(y * (num_bins_y - 1));
        if(x_bin < 0 || x_bin >= num_bins_x){continue;} 
        if(y_bin < 0 || y_bin >= num_bins_y){continue;} 
        histo[y_bin][x_bin] += data;
        counts[y_bin][x_bin] += 1;
    }

    for(int i = 0; i < num_bins_y; i++){
        for(int j = 0; j < num_bins_x; j++){
            if(i==0 && j==0){
            min_data = 1000000; 
            max_data = histo[0][0];
            }
            else{
            if(histo[i][j] < min_data && histo[i][j] != 0){min_data = histo[i][j];}
            if(histo[i][j] > max_data){max_data = histo[i][j];}
            } 
        }
    }
    cout << "max data " << max_data << " " << min_data << endl; 
    int num_pix_x = num_bins_x; 
    int num_pix_y = num_bins_y;
    unsigned char **RGB_data = new unsigned char*[num_bins_y]; 
    for(int i = 0; i < num_bins_y; i++)
        RGB_data[i] = new unsigned char[3*num_bins_x];
    cout << "min/max " << min_data << " " << max_data << endl;
    cout << "log min/max" << log10(min_data) << " " << log10(max_data) << endl;
    for(int i = 0; i < num_pix_y; i++){
        for(int j = 0; j < 3*num_pix_x; j++){
            RGB_data[i][j] = 0;
        }
    }
    double data_pts[5];
    if(log == 1){
    for(int ind = 0; ind < 5; ind++){
        data_pts[ind] = log10(min_data) + ((float) ind /  4.0) * (log10(max_data) - log10(min_data));
    }
    cout << data_pts[2] << endl;
    }
    else{
    for(int ind = 0; ind < 5; ind++){
        data_pts[ind] = min_data + ((float) ind /  4.0) * (max_data - min_data);
    }

    }
    unsigned int RGB_vals[5][3] = {{0,0,0}, {0,0,100}, {0,170, 170}, {170, 170, 170}, {255,255,255}};
    for(int i = 0; i < num_bins_y; i++){
        for(int j = 0; j < num_bins_x; j++){
            if(0){ //counts[i][j] for interpolation
                int close_i = 0; 
                int close_j = 0; 
                for(int l = 1; l < 1000; l++){
                    int min_val = -1*l; 
                    int max_val = l; 
                    int i_found[4*l], j_found[4*l], dist[4*l];
                    int found = 0;
                    for(int t = -1*l; t <= l; t++){
                        if(counts[min_val][t] != 0 && min_val > 0 && min_val < num_bins_x && t > 0 && t < num_bins_y){
                            i_found[found] = min_val; 
                            j_found[found] = t;
                            dist[found] = min_val*min_val + t*t; 
                            found++; 
                        }
                        if(counts[max_val][t] != 0 && max_val > 0 && max_val < num_bins_x && t > 0 && t < num_bins_y){
                            i_found[found] = max_val; 
                            j_found[found] = t;
                            dist[found] = max_val*max_val + t*t; 
                            found++; 
                        }
                    }
                    for(int t = 1-l; t <= l-1; l++){
                        if(counts[t][min_val] != 0 && t > 0 && t < num_bins_x && min_val > 0 && min_val < num_bins_y){
                            i_found[found] = t; 
                            j_found[found] = min_val;
                            dist[found] = min_val*min_val + t*t; 
                            found++; 
                        }
                    }
                    for(int t = 1-l; t <= l-1; l++){
                        if(counts[t][max_val] != 0 && t > 0 && t < num_bins_x && max_val > 0 && max_val < num_bins_y){
                            i_found[found] = t; 
                            j_found[found] = max_val;
                            dist[found] = max_val*max_val + t*t; 
                            found++; 
                        }
                    }
                    if(found > 0){
                        if(found > 1){
                            float min_dist = dist[0];
                            int close_i = i_found[0]; 
                            int close_j = j_found[0];
                            for(int c = 1; c < found; c++){
                                if(dist[c] < min_dist){
                                    min_dist = dist[c];
                                    close_i = i_found[c]; 
                                    close_j = j_found[c];
                                }
                            }
                        }
                        break;
                    }
                }
            int img_x_ind = (static_cast<float>(j) / static_cast<float>(num_bins_x)) * num_pix_x; 
            int img_y_ind = (static_cast<float>(i) / static_cast<float>(num_bins_y)) * num_pix_y;
            int lower_ind = 0; 
            int upper_ind = 1;
            for(int k = 1; k < 5; k++){
                if(log10(histo[close_i][close_j]) >= data_pts[lower_ind] && log10(histo[close_i][close_j])   < data_pts[upper_ind]){break;}
                else{
                    lower_ind++; 
                    upper_ind++; 
                }
            }
            float delta_l = (histo[close_i][close_j] - data_pts[lower_ind]);
            if(histo[close_i][close_j] == 0){delta_l = 0;}
            float RGB_val[3];
            float delta = 0.25; 
            for(int k = 0; k < 3; k++){
            RGB_val[k] = (RGB_vals[upper_ind][k] - RGB_vals[lower_ind][k])/(delta) *delta_l + RGB_vals[lower_ind][k];
            RGB_data[img_y_ind][3*(img_x_ind) + k] = static_cast<unsigned char>(RGB_val[k]);//changed from RGB_vak[k] to counts for bound check
            }
            counts[i][j] += 1;
            histo[i][j] = histo[close_i][close_j];
            }
            else{
            //int img_x_ind = (static_cast<double>(j) / static_cast<double>(num_bins_x)) * num_pix_x;
            //int img_y_ind = (static_cast<double>(i) / static_cast<double>(num_bins_y)) * num_pix_y;
            int img_x_ind = j; 
            int img_y_ind = i; 
            //I need to nearest neighbor fill pixels that have 0 count in to make full image with less refined areas
            int lower_ind = 0; 
            int upper_ind = 1;
            double delta_l = 0;
            if(log == 1){
            for(int k = 1; k < 4; k++){
                if(log10(histo[i][j]) >= data_pts[lower_ind] && log10(histo[i][j]) < data_pts[upper_ind]){break;}
                else{
                    lower_ind++; 
                    upper_ind++; 
                }
            }
            delta_l = (log10(histo[i][j]) - data_pts[lower_ind]);
            if(histo[i][j] == 0){delta_l = 0; lower_ind = 0;}
            }
            else{
            for(int k = 1; k < 4; k++){
                if(histo[i][j] >= data_pts[lower_ind] && histo[i][j] < data_pts[upper_ind]){break;}
                else{
                    lower_ind++; 
                    upper_ind++; 
                }
            }
            delta_l = (histo[i][j] - data_pts[lower_ind]);
            if(histo[i][j] == 0){delta_l = 0; lower_ind = 0;}
            }
            double RGB_val[3];
            double delta = data_pts[upper_ind] - data_pts[lower_ind];
            for(int k = 0; k < 3; k++){
            RGB_val[k] = (RGB_vals[upper_ind][k] - RGB_vals[lower_ind][k])/(delta) *delta_l + RGB_vals[lower_ind][k];
            RGB_data[img_y_ind][3*(img_x_ind) + k] = static_cast<unsigned char>(RGB_val[k]); //RGB_val[k] changed x and y image inds
            }
            }
        }
    }
    cout << "past setting rgb vals" << endl;
    int height = num_pix_y; 
    int width = num_pix_x; 
    //Now make the png file 
    FILE* fp = fopen("test.png", "wb"); 
    if(fp == NULL){
        cout << "Error making .png file" << endl;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); 
    if(png_ptr == NULL){
        cout << "Error with PNG ptr" << endl;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr); 
    if(info_ptr == NULL){
        cout << "Error with png info ptr" << endl;
    }
    png_init_io(png_ptr, fp); 
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT); 
    png_write_info(png_ptr, info_ptr); 
    for(int i = 0; i < height; i++){
        png_write_row(png_ptr, (png_bytep)RGB_data[i]);
    }
    cout << "done writing rows" << endl;
    png_write_end(png_ptr, NULL); 
    cout << "wrote end" << endl;
    fclose(fp);
    cout << "closed file ptr" << endl;

    for(int i = 0; i < num_bins_y; i++){
        delete[] histo[i]; 
        delete[] counts[i];
        delete[] RGB_data[i];
    }
    delete[] histo; 
    delete[] counts;
    delete[] RGB_data;
    return 1; 
}

int GetNumNotRefinedGrids(grid *grids, int num_grids){
    int numNotRefined = 0; 
    for(int i = 0; i < num_grids; i++){
        if((grids + i)->GetNumNotRefined() == 0 && (grids+i)->GetGridID()>0){
            numNotRefined++; 
        }
    }
    return numNotRefined; 
}

