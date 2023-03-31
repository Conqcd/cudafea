//  solver class

#include "vFESolver.hpp"
#include "Math/Solver.hpp"
#include<chrono>
#include<thread>
#include<cmath>

//============ Constructor / Destructor

// COLLECTIVE
vFESolver::vFESolver() 
{
    // set all flags to false
    Init();
}


vFESolver::~vFESolver()
{}

//============ Initialise

// initialise variables here
void vFESolver::Init()
{
    // setup and rebuild flags
    VOXELSIZE_DONE = false;
    VOXELDIMS_DONE = false;
    MATERIALS_DONE = false;
    MODEL_DONE = false;
    CONSTRAINTS_DONE = false;
    
    SOLVE_DONE = false;
    FINISH_DONE = false;
    
    // set default script/model file format
    SCRIPTVERSION = 2;
    
    
}// Init()

// check all minimal setup commands have been called. Basic sanity check before solve called
bool vFESolver::IsSetupDone()
{
    bool OK = VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE && MODEL_DONE && CONSTRAINTS_DONE;
    return OK;
}// IsSetupDone()

//========= Setters

// Sets the voxel dimensions and scale factor
bool vFESolver::SetVoxelSize(const double a, const double b, const double c, const double sf )
{
    bool OK(false);
    if(a * b * c * sf)
    { // check if not zero
        
        A_SIZE = a;
        B_SIZE = b;
        C_SIZE = c;
        SCALE_FACTOR = sf;
        VOXELSIZE_DONE = true;
        OK = true;
    }
    else
    {
        printf("ERROR: could not add voxelsizes and/or scale factor \n");
        // OK = false;
    }
    return OK;
}// SetVoxelSize()


// Sets the model dimensions in voxels : no default
bool vFESolver::SetVoxelDimensions(const xyzType dimx, const xyzType dimy, const xyzType dimz)
{
    bool OK(false);
    if(dimx * dimy * dimz)
    {
        DX = dimx;
        DY = dimy;
        DZ = dimz;
        DTOT = DX*DY*DZ;
        DXY  = DX*DY;
        NX   = DX + 1;
        NY   = DY + 1;
        NZ   = DZ + 1;
        NTOT = NX*NY*NZ;
        NXY  = NX*NY;
        OK   = true;
        VOXELDIMS_DONE = true;
        OK = true;
    }
    else
    {
        printf("ERROR: voxel dimensions must be non-zero\n");
        // OK = false;
    }
    return  OK;
}// SetVoxelDimensions()


// Sets the model dimensions in voxels : no default
bool vFESolver::SetTotElems(const xyzType elemtot)
{
    bool OK(false);
    if( elemtot )
    {
        ELEMTOT = elemtot;
        OK = true;
    }
    else
    {
        printf("ERROR: total number of elements must be non-zero\n");
        // OK = false;
    }
    return  OK;
}// SetTotElems()



bool vFESolver::SetMaxIter(const long int themaxiter)
{
    bool OK(false);
    MAXITER = themaxiter;
    return (OK = true);
}// SetMaxIter()

// set tolerance
bool vFESolver::SetTolerance(const double thetolerance)
{
    bool OK(false);
    TOLERANCE = thetolerance;
    return (OK = true);
}// SetTolerance()

bool vFESolver::SetAlgorithmFEA(const char *thesolvertype, const char *thepctype)
{
    bool OK(false);
    strcpy(SOLVERTYPE, thesolvertype);
    strcpy(PCTYPE, thepctype);
    return (OK = true);
}

//================ Loaders

// Create Material map
bool vFESolver::LoadMaterials(const char *filename)
{
    strcpy( MATERIALSFILE, filename );
    printf("material file is %s \n", filename);
    
    FILE * materialFile;
    
    midxType material_index;
    double youngs_mod;
    double poissons_rat;
    int remodel_flag;
    bool OK(true);
    
    if(VOXELSIZE_DONE)
    {
        printf("Creating Material Map \n");
        
        // compute gradient matrices : needed for computing LSM of each material
        for(int n = 0; n < SAMPLES_PER_ELEMENT; ++n){
    
            gradientMtx[n] = GradientMatrix(A_SIZE, B_SIZE, C_SIZE, n);
        }
        
        printf("Opening material file \n");
        materialFile = fopen(filename, "r");
        errno = 0;
        if(materialFile == NULL)
        {
            printf("ERROR: could not open material file \n");
            OK = false;
            return OK;
        }else
        {
            printf("material file open \n");
            
            while(fscanf(materialFile, "%d %lf %lf %d \n", &material_index, &youngs_mod, &poissons_rat, &remodel_flag) != EOF && OK)
            {// while
                
                OK = AddMaterial(material_index, youngs_mod, poissons_rat, remodel_flag);
                
            }// while adding materials
            if(!OK){
                printf("ERROR: couldn't add material \n");
                fclose(materialFile);
                OK = false;
                return OK;
            }// if material not added
            
            fclose(materialFile);
            
            MATERIALS_DONE = true;
        }
        
    }// if voxel sizes specified
    else
    {
        printf("ERROR: voxel sizes must be specified before loading materials \n");
        OK = false;
    }// voxel size or voxel dimensions not specified
    const int thread_num = std::thread::hardware_concurrency();
    const int max_threads = std::min((int)MateM.size(),thread_num);
    std::vector<std::thread> threads(max_threads);
    for (int i = 0; i < max_threads; i++)
    {
        threads[i] = std::thread([&,i](){
            for (int j = i; j < MateM.size(); j += max_threads)
            {
                MateM[j].lsm.create(gradientMtx);
            }
        });
    }
    for(auto& entry:threads)
    {
        entry.join();
    }
    return OK;
}// LoadMaterials()


// Load elements, create element and node sets
bool vFESolver::LoadModel(const char *filename)
{
    strcpy(MODELFILE, filename);

    midxType                    material;
    xyzType                     element_num, x, y, z;
    FILE *                      modelFile;
    bool OK(true);
    
    OK = VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE;
    if(!OK)
    {
        printf("ERROR: voxel size, dimensions and materials must be specified before loading model \n");
        return OK;
    }
    
    printf("Creating Element and Node Sets\n");
    
    modelFile = fopen(filename,"r");
    
    if(modelFile == NULL)
    {
        printf ("ERROR: could not open model file \n");
        OK = false;
        return OK;
    }
    else
    {
        // if default script/model file formats
        // read first two lines and continue
        if( SCRIPTVERSION != 2)
        {
         fscanf(modelFile, " ");
         fscanf(modelFile, "%d", &element_num);
            if(element_num)
                ELEMTOT = element_num;
            else
            {
                OK = false;
                printf ("ERROR: total number of elements must be non-zero \n");
                return OK;
            }
        }// if default script/model file formats
        
        while(fscanf(modelFile, "%d %d %d %d %d ", &element_num, &material, &x, &y, &z)  != EOF && OK)
        {
            // script coords are 1-based, solver coords are 0-based
            //--x; --y; --z;
            
			OK = AddElement(material, x, y, z);
        }
        if(!OK)
        {
            return OK;
        }
    }// else if file exists
    fclose (modelFile);
    
    OK = RenumberNodes();
	
    if(OK)
    {
        MODEL_DONE = true;
    }
    return OK;
}// LoadModel()


bool vFESolver::ComsolResult2Constraint(const char *filename)
{
    std::ifstream f;
    char buff[256];
    VecString lines;
    VecString_it vit = lines.begin();
    printf("Transfer File\n");
    try
    {
        f.open(filename);
        size_t cstart;
        std::string line;
        getline(f, line);getline(f, line);getline(f, line);
        getline(f, line);getline(f, line);getline(f, line);
        getline(f, line);getline(f, line);getline(f, line);
        for(std::string line; getline(f, line);)
        {
            cstart = line.find_first_not_of(" \t\n\r");
            
            if( cstart != std::string::npos )  lines.push_back(line);
        }
    }// try
    catch(std::exception &e)
    {
        printf("ERROR: could not open Comsol result file \n");
        return false;
    }// catch
    f.close();
    double x, y, z;
    double vx, vy, vz;
    std::ofstream consfile;
    consfile.open("temp/constraint2.txt");
    if(!consfile.is_open())
    {
        printf("ERROR: could not write output constraint file \n");
        return false;
    }
    for(auto& line:lines)
    {
        sscanf(line.c_str(), "%lf %lf %lf %lf %lf %lf", &x, &y, &z, &vx, &vy, &vz);
        if(vx == 0 && vy == 0 && vz == 0)
        {
            sprintf(buff,"SELECT_NODE_3D %d %d %d %lf %lf %lf 1 1 1\n",(int)x,(int)y,(int)z,vx,vy,vz);
            Node<xyzType> * tmpnode = new Node<xyzType>();
    
    
            tmpnode->x = x; tmpnode->y = y; tmpnode->z = z;

            auto nitr = NodeS.find(tmpnode);
            if(nitr == NodeS.end())
            {
                consfile.close();
                return false;
            }
            consfile << buff;
        }
    }
    consfile.close();
    return true;
}

// Read in constraint and force data
bool vFESolver::LoadConstraints(const char *filename){
    
    strcpy(CONSTRAINTSFILE, filename);
    
    bool OK(true);
    xyzType x, y, z;
    int cx, cy, cz;
    double vx, vy, vz;
    char str[80];
    
    VecString lines;
    VecString_it vit = lines.begin();
    NodeSet_it nodeptr = NodeS.begin();
    
    idxType totlines;
    idxType conscount(0);
    
    
    OK = ( VOXELSIZE_DONE && VOXELDIMS_DONE && MATERIALS_DONE && MODEL_DONE );
    if(!OK)
    {
        printf("ERROR: model parameters, materials and model need to be specified before constraints \n");
        return OK;
    }// if voxel sizes etc. specified
    
    
    std::ifstream consFile(filename);
    printf("Reading constraints\n");
    try
    {
        size_t cstart;
        for(std::string line; getline(consFile, line);)
        {
            // check if empty line
            cstart = line.find_first_not_of(" \t\n\r");
            
            // only add line if non-white space characters found
            if(cstart != std::string::npos)  lines.push_back(line);
        }
    }// try
    catch(std::exception &e)
    {
        printf("ERROR: could not open constraint file \n");
        OK = false;
        return OK;
    }// catch
    consFile.close();
    totlines = lines.size();

    
    Node<xyzType> * tmpnode = new Node<xyzType>();
	
    totalrhs = 0;
	
    for(vit = lines.begin(); vit != lines.end(); ++vit)
    {
        sscanf((*vit).c_str(), "%s %u %u %u %lf %lf %lf %d %d %d",&str, &x, &y, &z, &vx, &vy, &vz, &cx, &cy, &cz);

        totalrhs += !cx;
        totalrhs += !cy;
        totalrhs += !cz;
        //if(cx || cy || cz)
        //{
        //    //cx = !cx;  cy = !cy;  cz = !cz;
        //    totalrhs += !cx;
        //    totalrhs += !cy;
        //    totalrhs += !cz;
        //    // add constraint
        //    ++conscount;
        //    if (!strcmp(str, "SELECT_NODE_3D"))
        //        OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 0);
        //    else if (!strcmp(str, "PRESERVE_NODE_3D"))
        //        OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 1);
        //    else
        //    {
        //        printf("%d: ERROR: constraint command not recognised : %s \n", MPIrank, str);
        //        OK = false;
        //        return OK;
        //    }
        //    if(!OK)
        //    {
        //        printf("%d: ERROR: could not add constrained node %d : [%d, %d, %d] \n", MPIrank, conscount, x, y, z);
        //        return OK;
        //    }
        //}
        //else
        //    totalrhs += 3;
    }// for line in lines

    totalrhs = DOF_3D * (totlines - conscount);
    // forcevalue.resize(totalrhs);
    // forcecolidx.resize(totalrhs);

	
    conscount = 0; // not important here
    for(vit = lines.begin(); vit != lines.end(); ++vit){// is a force
        sscanf ((*vit).c_str(), "%s %d %d %d %lf %lf %lf %d %d %d",&str, &x, &y, &z, &vx, &vy, &vz, &cx, &cy, &cz);
  
        //--x; --y; --z;
        //cx = !cx;  cy = !cy;  cz = !cz;
        
        // add force;
        ++conscount;
        
        if (!strcmp(str, "SELECT_NODE_3D"))
            OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 0);
        else if (!strcmp(str, "PRESERVE_NODE_3D"))
            OK = AddConstraint(x, y, z, vx, vy, vz, cx, cy, cz, 1);
        else
        {
            printf("ERROR: constraint command not recognised : %s \n", str);
            OK = false;
            return OK;
        }
        if(!OK)
        {
            printf("ERROR: could not add constrained node %d : [%d, %d, %d] \n", conscount, x, y, z);
            // return OK;
        }
        
    }// for lines

    OK = true;
    Idx_bias.resize(NodeS.size()* DOF_3D - totlines * DOF_3D + totalrhs);
	
    idxType num = 0;
    int csum = 0;
    NodeSet_const_it cnitr;

    for (cnitr = NodeS.begin(); cnitr != NodeS.end(); ++cnitr)
    {
    	if(!(*cnitr)->constraint)
        {
            Idx_bias[num - csum] = num;
            num++;
            Idx_bias[num - csum] = num;
            num++;
            Idx_bias[num - csum] = num;
            num++;
    		continue;
    	}
        if (!(*cnitr)->constraint->cx)
        {
            csum++;
        }else
        {
            Idx_bias[num - csum] = num;
        }
        num++;
        if (!(*cnitr)->constraint->cy)
        {
            csum++;
        }
        else
        {
            Idx_bias[num - csum] = num;
        }
        num++;
        if (!(*cnitr)->constraint->cz)
        {
            csum++;
        }
        else
        {
            Idx_bias[num - csum] = num;
        }
        num++;
    }
	
    CONSTRAINTS_DONE = true;
    
    return OK;
}// LoadConstraints()

//==================== Output

bool vFESolver::toNASTRAN(const char* filename)
{
    std::ofstream f;
    f.open(filename);
    if (!f.is_open())
        return false;
    f << "$ Generated by COMSOL 5.6.0.341\n$ Mar 2 2023, 14:28\nBEGIN BULK\n$ Grid data section\n";
    const char* grid = "GRID    ";
    char space[16] = "        ";
    for (auto& node : NodeS)
    {
        std::string line = grid;
        sprintf(space,"%8llu",node->idx + 1);
        line += space;
        sprintf(space,"        ");
        line += space;
        sprintf(space,"%8.6lf",node->x * A_SIZE * SCALE_FACTOR);
        space[8] = '\0';
        line += space;
        sprintf(space,"%8.6lf",node->y * B_SIZE * SCALE_FACTOR);
        space[8] = '\0';
        line += space;
        sprintf(space,"%8.6lf",node->z * C_SIZE * SCALE_FACTOR);
        space[8] = '\0';
        line += space;
        f << line << std::endl;
    }
    f << "$ Element data section\n";
    const char* chexa = "CHEXA   ";
    const char* cont1 = "+CONT";
    for (auto& ele : ElemS)
    {
        std::string line = chexa;
        sprintf(space,"%8llu",ele->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->material + 1);
        line += space;

        // ************* Node ****************
        sprintf(space,"%8u",ele->nodes[0]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[1]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[3]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[2]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[4]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[5]->idx + 1);
        line += space;
        sprintf(space,"%s\n",cont1);
        line += space;
        sprintf(space,"%s   ",cont1);
        line += space;
        sprintf(space,"%8u",ele->nodes[7]->idx + 1);
        line += space;
        sprintf(space,"%8u",ele->nodes[6]->idx + 1);
        line += space;
        f << line << std::endl;
    }
 
    f << "$ Property data section\n";

    const char* psolid = "PSOLID  ";

    for (auto& mat : MateM)
    {
        std::string line = psolid;
        sprintf(space,"%8llu",mat.idx + 1);
        line += space;
        f << line << std::endl;
    }
    space[8] = '\0';
    f << "ENDDATA";
    f.close();
    return true;
}

bool vFESolver::PrintDisplacements(const char *filename)
{
    // get displacements from node and output 
    
    strcpy(OUTPUTFILE, filename);
    
    bool OK(false);
    idxType globalID;
    
    xyzType nx, ny, nz;
    Scalar dx, dy, dz;
    FILE * outFile;
    NodeSet_it nitr;
    
    try
    {
        outFile = fopen(filename,"w");
        if(outFile != NULL){
            
            printf("output file is open for displacements\n");
        }else{
            printf("ERROR: output file NULL\n");
            return OK;
        }
        
        char materialFilename[] = "./axial_material.txt";
        fprintf(outFile,"nNodes: %d\n", NodeS.size());
        fprintf(outFile,"Nodes: %d  %d  %d\n", NX, NY, NZ);
        fprintf(outFile,"Materials: %s\n", MATERIALSFILE);
        
        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            nx = (*nitr)->x;
            ny = (*nitr)->y;
            nz = (*nitr)->z;
            dx = (*nitr)->dx;
            dy = (*nitr)->dy;
            dz = (*nitr)->dz;
            globalID = (*nitr)->idx;
			//globalID = nx + ny*NX + nz*NXY;      
      
			fprintf(outFile, "%d   %d   %d   %d   %.9le   %.9le   %.9le  \n", globalID, nx, ny, nz, dx, dy, dz);
        }// for each node in Node Set
        fclose(outFile);
        DISPLACEMENT_DONE = true;
        OK = true;
    }
    catch(std::exception &e)
    {
        printf("ERROR: could not print displacements \n");
        fclose(outFile);
        OK = false;
        DISPLACEMENT_DONE = false;
        return OK;
    }
    return (OK);
    
} // PrintDisplacements()

bool vFESolver::PrintStrain(const char *filename)
{
    // get strains from node and output 
    
    strcpy(OUTPUTFILE, filename);
    
    bool OK(false);
    idxType globalID;
    
    xyzType nx, ny, nz;
    Scalar dx, dy, dz;
    FILE * outFile;
    NodeSet_it nitr;

    ComputeStrainAndStrss();

    try
    {
        outFile = fopen(filename,"w");
        if(outFile != NULL){
            
            printf("output file is open for strains\n");
        }else{
            printf("ERROR: output file NULL\n");
            return OK;
        }
        
        fprintf(outFile,"nNodes: %d\n", NodeS.size());
        fprintf(outFile,"Nodes: %d  %d  %d\n", NX, NY, NZ);
        fprintf(outFile,"Materials: %s\n", MATERIALSFILE);
        
        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            nx = (*nitr)->x;
            ny = (*nitr)->y;
            nz = (*nitr)->z;
            dx = (*nitr)->dx;
            dy = (*nitr)->dy;
            dz = (*nitr)->dz;
            globalID = (*nitr)->idx;
      
			fprintf(outFile, "%d   %d   %d   %d   %.9le   %.9le   %.9le  \n", globalID, nx, ny, nz, dx, dy, dz);
        }// for each node in Node Set
        fclose(outFile);

        OK = true;
    }
    catch(std::exception &e)
    {
        printf("ERROR: could not print strains \n");
        fclose(outFile);
        OK = false;
        return OK;
    }
    return (OK);
}

bool vFESolver::PrintStress(const char *filename)
{
    // get strsses from node and output 
    
    strcpy(OUTPUTFILE, filename);
    
    bool OK(false);
    idxType globalID;
    
    xyzType nx, ny, nz;
    Scalar vm,fx,fy,fz,fxy,fzy,fzx;
    FILE * outFile;
    NodeSet_it nitr;

    try
    {
        outFile = fopen(filename,"w");
        if(outFile != NULL){
            
            printf("output file is open for strsses\n");
        }else{
            printf("ERROR: output file NULL\n");
            return OK;
        }
        
        fprintf(outFile,"nNodes: %d\n", NodeS.size());
        fprintf(outFile,"Nodes: %d  %d  %d\n", NX, NY, NZ);
        fprintf(outFile,"Materials: %s\n", MATERIALSFILE);
        
        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            nx = (*nitr)->x;
            ny = (*nitr)->y;
            nz = (*nitr)->z;
            vm = (*nitr)->von_mises / (*nitr)->ct;
            fx = (*nitr)->fx / (*nitr)->ct;
            fy = (*nitr)->fy / (*nitr)->ct;
            fz = (*nitr)->fz / (*nitr)->ct;
            fxy = (*nitr)->fxy / (*nitr)->ct;
            fzy = (*nitr)->fzy / (*nitr)->ct;
            fzx = (*nitr)->fzx / (*nitr)->ct;

            globalID = (*nitr)->idx;
      
			fprintf(outFile, "%d   %d   %d   %d   %.9le   %.9le   %.9le   %.9le   %.9le   %.9le   %.9le\n", globalID, nx, ny, nz, vm, fx, fy, fz, fxy, fzy, fzx);
        }// for each node in Node Set
        fclose(outFile);

        OK = true;
    }
    catch(std::exception &e)
    {
        printf("ERROR: could not print strsses \n");
        fclose(outFile);
        OK = false;
        return OK;
    }
    return (OK);
}

//==================== Adders

// Adds a material with given properties, computes and stores the LSM
bool vFESolver::AddMaterial(const midxType index, const double ym, const double pr, const int rf){
    bool OK(false);
    try
    {
        LocalStiffnessMatrix lsm(A_SIZE, B_SIZE, C_SIZE, ym, pr, gradientMtx);
         
        // Material * new_material = new Material(index, ym, pr, lsm, rf);
        
        if(MateM.size() <= index)
            MateM.resize(index + 1);
        MateM[index] = {index,ym,pr,lsm,(bool)rf};
        
        // delete new_material;
        // new_material = 0;
        
        OK = true;
    }
    catch (std::exception &e)
    {
        printf("ERROR: could not add material : %d %lf %lf %d\n", index, ym, pr, rf);
        OK = false;
    }
    return OK;
}


// Add element and associated nodes if they don't already exist
bool vFESolver::AddElement(const midxType material, const xyzType x, const xyzType y, const xyzType z){
    std::pair< ElementSet_it , bool>  elem_pair;
    std::pair< NodeSet_it    , bool>  node_pair;
    xyzType  nx, ny, nz;
    
    NodeSet_it                   nitr;
    ElementSet_it                eitr;
    
    idxType                      ecount(0), ncount(0);
    bool OK(true);
    
    // if material is in the MaterialMap i.e. an accepted material type
    try
    {
        // does material exist ?
        // if (!(MateM.find(material) != MateM.end()))
        if (material >= MateM.size())
        {
            printf("ERROR: material %d does not exist \n",material);
            OK = false;
            return OK;
        }

        Element<xyzType> * elem = new Element<xyzType>();
        elem->x = x; elem->y = y; elem->z = z;
        elem->material = material;
        elem->idx = x + y * DX + z * DXY;
        ++ecount;
        
        elem_pair = ElemS.insert(elem);
        eitr = elem_pair.first;
        
        Node<xyzType>* tmpnode = new Node<xyzType>();
        
        for (int n = 0; n < NODES_PER_ELEMENT; ++n) {// for each node
            
            //calculate nodal coordinates/idx from element coordinates
            nx = x + (n % 2);
            ny = y + (n / 2) * (n < 4) + (n / 6) * (n > 4);
            nz = z + (n / 4);
            tmpnode->x = nx;
            tmpnode->y = ny;
            tmpnode->z = nz;
            
            nitr = NodeS.find(tmpnode); // try to find node in Node Map
            
            if (nitr == NodeS.end()) { // if node not in Node Map

                Node<xyzType> * newnode = new Node<xyzType>();
                newnode->x = nx;
                newnode->y = ny;
                newnode->z = nz;
                newnode->idx = nx + ny * NX + nz * NXY;
    
                ++ncount;
                node_pair = NodeS.insert(newnode);
                nitr = node_pair.first;
                
            }// if node exists
            
            (*eitr)->nodes[n] = (*nitr);
            (*nitr)->elems[n] = (*eitr);
            
        }// for each node
        
        delete tmpnode;
        
    }// try
    catch(std::exception &e){
        printf("ERROR: could not add element \n");
        OK = false;
    }
    return OK;
}// AddElement()


// Adds constraint to the node
bool vFESolver::AddConstraint(const xyzType x, const xyzType y, const xyzType z, const double vx, const double vy, const double vz, const int cx, const int cy, const int cz, const int pf)
{
    
    NodeSet_it nitr = NodeS.begin();
    idxType nidx;
    Node<xyzType> * tmpnode = new Node<xyzType>();
    Constraint<xyzType> * cons = new Constraint<xyzType>();
    bool OK(true);
    
    cons->cx = !cx;
    cons->cy = !cy;
    cons->cz = !cz;
    cons->vx = vx;
    cons->vy = vy;
    cons->vz = vz;
    cons->preserve = pf;
    
    tmpnode->x = x; tmpnode->y = y; tmpnode->z = z;

    nitr = NodeS.find(tmpnode);

    
    if (nitr != NodeS.end())
    {
        cons->node = (*nitr);
        (*nitr)->constraint = cons;

        nidx = (*nitr)->idx;
    	
        if (!cx) // add force
        {
            // forcevalue[forcecount] = vx;
            forcevalue.push_back(vx);
            // forcecolidx[forcecount] = nidx * DOF_3D + 0;
            forcecolidx.push_back(nidx * DOF_3D + 0);
            ++forcecount;
        }
        if (!cy) // add force
        {
            // forcevalue[forcecount] = vy;
            forcevalue.push_back(vy);
            // forcecolidx[forcecount] = nidx * DOF_3D + 1;
            forcecolidx.push_back(nidx * DOF_3D + 1);
            ++forcecount;
        }
        if (!cz) // add force
        {
            // forcevalue[forcecount] = vz;
            forcevalue.push_back(vz);
            // forcecolidx[forcecount] = nidx * DOF_3D + 2;
            forcecolidx.push_back(nidx * DOF_3D + 2);
            ++forcecount;
        }
    	// if force
    }// if node exists
    else
    {
        delete cons;
        delete tmpnode;
        printf("ERROR: node R[%d, %d, %d], V[%f, %f, %f], C[%d, %d, %d] does not exists, cannot add constraint \n", x, y, z,vx,vy,vz,cx,cy,cz);
        OK = false;
        return OK;
    }
    return OK;
}// AddConstraint()


// Delete node
bool vFESolver::RemoveNode(const xyzType nx, const xyzType ny, const xyzType nz)
{
    
    bool OK(true);
    return OK;
}

// Delete an element to the model (along with its nodes)
bool vFESolver::RemoveElement(const midxType material, const xyzType x, const xyzType y, const xyzType z)
{
    bool OK(true);
    return OK;
}

// Delete constraint
bool vFESolver::RemoveConstraint(const xyzType x, const xyzType y, const xyzType z)
{
    
    bool OK(true);
    return OK;
}

//========== Getters

// Get a particular node by coords
Node<xyzType>* vFESolver::GetNode(const xyzType x, const xyzType y, const xyzType z)
{
    Node<xyzType>* tmpnode = new Node<xyzType>();
    tmpnode->x = x;
    tmpnode->y = y;
    tmpnode->z = z;
    NodeSet_it nodeptr =  NodeS.find(tmpnode);
    
    return (*nodeptr);
}// GetNode()

// Gets a particular element by coords
Element<xyzType>* vFESolver::GetElement(const xyzType x, const xyzType y, const xyzType z)
{
    return (*ElemS.begin());
}// GetElement()

// Gets a set of all the (local) elements
vFESolver::ElementSet vFESolver::GetLocalElements()
{
    return ElemS;
}// GetLocalElements()

// Gets a set of all the (local) nodes
vFESolver::NodeSet vFESolver::GetLocalNodes()
{
    return NodeS;
}// GetLocalNodes()


Constraint<xyzType>* vFESolver::GetConstraint(const xyzType x, const xyzType y, const xyzType z)
{
    Constraint<xyzType>* constraint = (GetNode(x, y, z))->constraint;
    return constraint;
}


//================== Allocate system matrix / vector rows

// allocate gsm rows and columns according to number of processes
bool vFESolver::AllocateLocalMatrix(Matrix& gsm)
{
    printf("NodeS.size() = %d \n", NodeS.size());
    
    
    // IF PARALLEL
    double gr(globalgsmrows);
    double gc(globalgsmcols);
    // localgsmrows = floor(gr / MPIsize);
    // localgsmcols = floor(gc / MPIsize);
    
    // if you are the top ranking process
    // take remaining rows/columns
    // if (MPIrank == MPIsize - 1)
    // {
    //     localgsmrows = gr - (MPIsize - 1) * floor(gr / MPIsize);
    //     localgsmcols = gc - (MPIsize - 1) * floor(gc / MPIsize);
    // }
    
    diagonalterms = NUM_TERMS;
    offdiagonalterms = NUM_TERMS;
    
    localrhslength = localgsmcols;
    localsollength = localgsmcols;
    
    printf("GSM Alloc: lrows=%d, lcols=%d \n", localgsmrows, localgsmrows);

    // MatSetSizes(*gsm, localgsmrows, localgsmcols, globalgsmrows, globalgsmcols);
    gsm.reset(globalgsmrows,globalgsmcols);
    // MatSetType(*gsm, MATMPIAIJ);
    // MatMPIAIJSetPreallocation(*gsm, diagonalterms, PETSC_NULL, offdiagonalterms, PETSC_NULL);
    gsm.PreAllocation(diagonalterms);
        
    
    // check rows
    int grows, gcols;
    
    grows = gsm.get_row();
    gcols = gsm.get_col();
    
    printf("MatGetSize of GSM : rows=%d, cols=%d\n",grows, gcols);
    
    return 0;
}

// Allocate RHS rows according to number of processes
bool vFESolver::AllocateLocalVec(Vector& vec)
{
    // if(MPIsize == 1)
    {// IF SERIAL
        // VecSetSizes(*vec, globalgsmrows, globalgsmrows);
        // VecSetType(*vec, VECSEQ);
        // vec.reset(globalgsmrows);
    }
    // else
    {// IF PARALLEL
        // VecSetSizes(*vec, localgsmcols, globalgsmrows);
        // VecSetType(*vec, VECMPI);
        vec.reset(globalgsmrows);
    }
    return 0;
}


//=============== System matrices

// Build Global Stiffness Matrix
bool vFESolver::ComputeGSM(Matrix& GSM)
{
    printf("In GSM\n");
    
    globalgsmrows =  NodeS.size() * DOF_3D;
    globalgsmcols =  NodeS.size() * DOF_3D; // NUM_TERMS = max number of non-zero elements per row
    
    printf("GSM ROWS [COLS] = %d \n", globalgsmrows);
    
    // MatInfo matinfo; // use to get matrix info
    
    // allocate GSM rows
    vFESolver::AllocateLocalMatrix(GSM);
    
    // inclusive of first, exclusive of last

    const unsigned int lsmlen = NODES_PER_ELEMENT * DOF_3D; // 24
    const int thread_num = std::thread::hardware_concurrency();
    // const int thread_num = 1;
    const int max_threads = std::min((int)NodeS.size(),thread_num);

    
    auto startT = std::chrono::high_resolution_clock::now();
    // Threads
    std::vector<std::thread> threads(max_threads);
    for (int i = 0; i < max_threads; i++)
    {
        threads[i] = std::thread([&,i](){

        std::map<idxType, idxType> tmp_gsmcolidx;
        idxType currentcol(0);
        std::vector<idxType> gsmCol(NUM_TERMS);
        std::vector<idxType> gsmRowX(1), gsmRowY(1), gsmRowZ(1);
        
        std::vector<Scalar> gsmvalRowX(NUM_TERMS); // row-centric storage etc.
        std::vector<Scalar> gsmvalRowY(NUM_TERMS);
        std::vector<Scalar> gsmvalRowZ(NUM_TERMS);
        
        bool consnodeX(1); // is curnode constrained in X dim Default is 1 (true)
        bool consnodeY(1); // is curnode constrained in Y dim Default is 1 (true)
        bool consnodeZ(1); // is curnode constrained in Z dim Default is 1 (true)
        
        bool consnodeNeighbourX(1); // is curennode constrained in X dim Default is 1  (true)
        bool consnodeNeighbourY(1); // is curennode constrained in X dim Default is 1  (true)
        bool consnodeNeighbourZ(1); // is curennode constrained in X dim Default is 1  (true)
        
        bool constraintX[DOF_3D] = {1,1,1}; // should X dim term be added to GSM Default is 1 (yes)
        bool constraintY[DOF_3D] = {1,1,1}; // should Y dim term be added to GSM Default is 1 (yes)
        bool constraintZ[DOF_3D] = {1,1,1}; // should Z dim term be added to GSM Default is 1 (yes)
        
        ////////////////////////
        // BUILD GSM BY LOOPING THROUGH NODES IN NodeS
        ///////////////////////
        
        idxType renumNode(0), renumNodeNeighbour(0), nodeL(0), nodeNeighbourL(0);
        idxType e, enn;
        midxType materialidx;
        xyzType c, colindex;
        auto cnitr = NodeS.begin();
        for (int j = 0; j < i; j++)
        {
            ++cnitr;
        }
        
        for(;cnitr != NodeS.end();)
        {
            renumNode = vFESolver::GetNodeIndex(cnitr); // renumbered or local index
            
            // reset stuff...
            tmp_gsmcolidx.clear();
            
            // for (c = 0; c < NUM_TERMS; ++c) {
            //     gsmCol[c] = 0;
            //     gsmvalRowX[c] = 0;
            //     gsmvalRowY[c] = 0;
            //     gsmvalRowZ[c] = 0;
            // }
            gsmCol.clear();
            gsmvalRowX.clear();
            gsmvalRowY.clear();
            gsmvalRowZ.clear();
            
            gsmRowX[0] = renumNode * DOF_3D  + 0;
            gsmRowY[0] = renumNode * DOF_3D  + 1;
            gsmRowZ[0] = renumNode * DOF_3D  + 2;
            
            Constraint<xyzType>* nodecons = vFESolver::GetNodeCons(cnitr);
            if(nodecons){
                consnodeX = nodecons->cx;
                consnodeY = nodecons->cy;
                consnodeZ = nodecons->cz;
            } else{
                consnodeX = 1;
                consnodeY = 1;
                consnodeZ = 1;
            }
            
            // enn = neighbour, e = elem
            for(int elem = 0; elem < ELEMENTS_PER_NODE; ++elem){// for each element the node belongs to (max 8)
                
                if(vFESolver::GetNodeElement(cnitr, elem)){
                    materialidx = vFESolver::GetElementMaterial(cnitr, elem);
                    nodeL = elem;
                    
                    for(int neighbour = 0; neighbour < NODES_PER_ELEMENT; ++neighbour){// for each neighbouring node on element e
                        renumNodeNeighbour = vFESolver::GetNodeNeighbourIndex(cnitr, elem, neighbour); //get renumbered index
                        
                        if(tmp_gsmcolidx.find(renumNodeNeighbour) == tmp_gsmcolidx.end()){// if neighbouring doesn't already have a column number
                            tmp_gsmcolidx[renumNodeNeighbour] = gsmCol.size() / 3;
                            gsmCol.push_back(0);
                            gsmCol.push_back(0);
                            gsmCol.push_back(0);
                            gsmvalRowX.push_back(0);
                            gsmvalRowX.push_back(0);
                            gsmvalRowX.push_back(0);
                            gsmvalRowY.push_back(0);
                            gsmvalRowY.push_back(0);
                            gsmvalRowY.push_back(0);
                            gsmvalRowZ.push_back(0);
                            gsmvalRowZ.push_back(0);
                            gsmvalRowZ.push_back(0);
                            // tmp_gsmcolidx[renumNodeNeighbour] = gsmcolcount;
                            // ++gsmcolcount;
                        }
                        currentcol = tmp_gsmcolidx[renumNodeNeighbour];
                        
                        Constraint<xyzType>* neighbourcons = vFESolver::GetNodeNeighbourCons(cnitr, elem, neighbour);
                        if(neighbourcons){
                            consnodeNeighbourX = neighbourcons->cx; // constrained //
                            consnodeNeighbourY = neighbourcons->cy; // constrained //
                            consnodeNeighbourZ = neighbourcons->cz; // constrained //
                        }else{
                            consnodeNeighbourX = 1;
                            consnodeNeighbourY = 1;
                            consnodeNeighbourZ = 1;
                        }
                        
                        nodeNeighbourL = neighbour; // local index of neighbouring node on element e
                        
                        constraintX[0] = consnodeX && consnodeNeighbourX;
                        constraintX[1] = consnodeX && consnodeNeighbourY;
                        constraintX[2] = consnodeX && consnodeNeighbourZ;
                        
                        constraintY[0] = consnodeY && consnodeNeighbourX;
                        constraintY[1] = consnodeY && consnodeNeighbourY;
                        constraintY[2] = consnodeY && consnodeNeighbourZ;
                        
                        constraintZ[0] = consnodeZ && consnodeNeighbourX;
                        constraintZ[1] = consnodeZ && consnodeNeighbourY;
                        constraintZ[2] = consnodeZ && consnodeNeighbourZ;
                        
                        if(nodeL == nodeNeighbourL){
                            if(!consnodeX){
                                constraintX[0] = 1; constraintX[1] = 0; constraintX[2] = 0;
                            }// if consnodeX
                            if(!consnodeY){
                                constraintY[0] = 0; constraintY[1] = 1; constraintY[2] = 0;
                            }// if consnodeY
                            if(!consnodeZ){
                                constraintZ[0] = 0; constraintZ[1] = 0; constraintZ[2] = 1;
                            }// if consnodeZ
                        }// if diagonal element of GSM
                        
                        for (c = 0; c < 3; ++c){// for each x,y,z component/column
                            colindex = currentcol * DOF_3D + c;
                            gsmCol[colindex] = renumNodeNeighbour * DOF_3D + c;
                            
                            int xrowidx = (nodeL * DOF_3D + 0) * lsmlen + nodeNeighbourL * DOF_3D + c;
                            int yrowidx = (nodeL * DOF_3D + 1) * lsmlen + nodeNeighbourL * DOF_3D + c;
                            int zrowidx = (nodeL * DOF_3D + 2) * lsmlen + nodeNeighbourL * DOF_3D + c;
                            
                            // add gsm value. Access correct LSM in MATERIALMAP using material index materialidx
                            gsmvalRowX[colindex] += vFESolver::GetLSMValue(materialidx, xrowidx) * constraintX[c];  // x row
                            gsmvalRowY[colindex] += vFESolver::GetLSMValue(materialidx, yrowidx) * constraintY[c];  // y row
                            gsmvalRowZ[colindex] += vFESolver::GetLSMValue(materialidx, zrowidx) * constraintZ[c];  // z row
                        }// for each component
                        
                    }// for each neighbouring node enn on element e
                }// if valid element e
            }// for each element e of current node
            
            GSM.insertValues(gsmRowX,gsmCol,gsmvalRowX);
            GSM.insertValues(gsmRowY,gsmCol,gsmvalRowY);
            GSM.insertValues(gsmRowZ,gsmCol,gsmvalRowZ);
            for(int j = 0; j < max_threads && cnitr != NodeS.end();cnitr++,j++);
        }// for each node
        });
    }
    
    // Get info about GSM allocation etc.
    // MatInfo info;
    double  mal = 0, nz_a = 0, nz_u = 0, nz_un = 0;
    
    // MatGetInfo(*GSM, MAT_LOCAL, &info);
    // mal  = info.mallocs;
    // nz_a = info.nz_allocated;
    // nz_u = info.nz_used;
    // nz_un = info.nz_unneeded;
    

    for(auto& entry:threads)
    {
        entry.join();
    }

    auto endT = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endT - startT);
    printf("GSM Building Time = %.5le ms\n", static_cast<double>(duration.count()));
    printf("GSM local info : mal = %lf, non-zero_allocated = %lf, non-zero_used = %lf, non-zero_unneeded = %lf \n", mal, nz_a, nz_u, nz_un);
    printf("Leaving GSM\n");
    return true;
}// computeGSM()



// Build RHS force vector
bool vFESolver::ComputeRHS(Vector& rhs){
    
    int xs,ys,zs,nx,ny,nz;
    int ii,jj,kk;
    int size;
    
    printf("In computeRHS\n");
    
    
    AllocateLocalVec(rhs);
    
    // VecGetSize(*rhs, &size);
    size = rhs.size();
    
    printf("VecGetSize : size = %d\n",size);
    
    // VecSet(*rhs,0);
    rhs.fill(0);
    // VecSetValues(*rhs, totalrhs, forcecolidx, forcevalue, INSERT_VALUES);
    rhs.setvalues(forcecolidx, forcevalue);
    
    printf("Leaving RHS \n");
    
    return true;
    
}// ComputeRHS()

void PrintMatrix(const Matrix& mat)
{
    for (int i = 0; i < mat.get_row(); i++)
    {
        for (int j = 0; j < mat.get_col(); j++)
        {
            std::cout << mat.index(i,j) << "\t";
        }
        std::cout << std::endl;       
    }
    
}

void ComputeB(std::vector<DenseMatrix>& Bs)
{
    
    const Scalar b[SAMPLES_PER_ELEMENT][6][24] = {
        -1	,0	,0	,1	,0	,0	,0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,-1,	0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0,	0,
        0	,0	,-1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0,	0,
        -1,	-1,	0	,0	,1	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0,	0,
        0	,-1	,-1	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0,	0,
        -1,	0	,-1	,0	,0	,1	,0	,0	,0	,0	,0	,0	,1	,0	,0	,0	,0	,0	,0	,0	,0	,0	,0,	0,
        
        -1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	0,	0,	0	,-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,
        0	,-1,	0,	-1,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	0,	0,	-1	,-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	-1,	-1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,

        0	,0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	0,	0,	0,	0,	0	,0	,-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,
        -1,	0,	0,	0,	0,	0,	1	,-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0	,0,	-1,	0,	0,	0,	0	,-1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,
        0	,0,	0,	0,	0,	0,	-1	,0	,-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,

        0,	0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,
        0,	0,	0,	-1,	0,	0,	0,	-1,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	-1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	-1,	-1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,

        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,
        0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	-1,	0,	0,	1,	0,	1,	0,	0,	0,	0,	0,
        0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,
        -1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,

        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,
        0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	-1,	1,	0,	0,	0,	0,	1,	0,	0,
        0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	-1,	0,	0,	0,	0,	0,	1,
        0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	1,	0,	1,	0,	0,	0,	0,	0,	0,

        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	-1,	0,	0,	1,	0,
        0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	1,	1,	0,	0,	0,
        0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	-1,	0,	0,	1,

        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	1,	0,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	1,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	-1,	0,	1,	1,	0,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	1,	1,
        0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-1,	1,	0,	1,
    };
    for (int i = 0; i < SAMPLES_PER_ELEMENT; i++)
    {
        for (int j = 0; j < Bs[i].get_row(); j++)
        {
            for (int k = 0; k < Bs[i].get_col(); k++)
            {
                Bs[i].insert(j,k,b[i][j][k]);
            }
        }
    }
    
}

bool vFESolver::ComputeStrainAndStrss()
{
    if(!DISPLACEMENT_DONE) return false;
    std::vector<DenseMatrix> Bs(NODES_PER_ELEMENT,{6,24});
    ComputeB(Bs);

    for (auto eitr = ElemS.begin(); eitr != ElemS.end(); ++eitr)
    {
        // Vector strain_gauss[SAMPLES_PER_ELEMENT],stress_gauss[SAMPLES_PER_ELEMENT];
        Vector strain[NODES_PER_ELEMENT],stress[NODES_PER_ELEMENT],mStrain[NODES_PER_ELEMENT],mStress[NODES_PER_ELEMENT];
        Vector d(NODES_PER_ELEMENT * 3);
        for (int i = 0; i < NODES_PER_ELEMENT; i++)
        {
            auto node = (*eitr)->nodes[i];
            d.setvalue(i * 3 + 0,node->dx);
            d.setvalue(i * 3 + 1,node->dy);
            d.setvalue(i * 3 + 2,node->dz);
        }
        Vector strain_Gavg(6),stress_Gavg(6);
        for (int i = 0; i < SAMPLES_PER_ELEMENT; i++)
        {
            // strain_gauss[i] = gradientMtx[i].getmat() * d;
            strain[i] = Bs[i] * d;
            // PrintMatrix(gradientMtx[i].getmat());
            // strain_Gavg += strain_gauss[i];
            // strain_gauss[i].scale(gradientMtx[i].getJdet());
            stress[i] = MateM[(*eitr)->material].lsm.propMtx.getmat() * strain[i];
            // PrintMatrix(MateM[(*eitr)->material].lsm.propMtx.getmat());
            // stress_Gavg += stress_gauss[i];
            // stress[i] = stress_gauss[i];
        }
        // strain_Gavg.scale(1.0 / SAMPLES_PER_ELEMENT);
        // stress_Gavg.scale(1.0 / SAMPLES_PER_ELEMENT);
        for (int i = 0; i < NODES_PER_ELEMENT; i++)
        {
            mStrain[i] = Math::solve3equation(1.0,-strain[i][0] - strain[i][1] - strain[i][2],
            (strain[i][0] * strain[i][1] + strain[i][2] * strain[i][1] + strain[i][0] * strain[i][2] 
            - strain[i][3] * strain[i][3] - strain[i][4] * strain[i][4] - strain[i][5] * strain[i][5]),
            -strain[i][0] * strain[i][1] * strain[i][2] - 2 * strain[i][3] * strain[i][4] * strain[i][5]
            + strain[i][0] * strain[i][4] * strain[i][4] + strain[i][1] * strain[i][5] * strain[i][5] +
            strain[i][2] * strain[i][3] * strain[i][3]);
            mStress[i] = Math::solve3equation(1.0,-stress[i][0] - stress[i][1] - stress[i][2],
            (stress[i][0] * stress[i][1] + stress[i][2] * stress[i][1] + stress[i][0] * stress[i][2] 
            - stress[i][3] * stress[i][3] - stress[i][4] * stress[i][4] - stress[i][5] * stress[i][5]),
            -stress[i][0] * stress[i][1] * stress[i][2] - 2 * stress[i][3] * stress[i][4] * stress[i][5]
            + stress[i][0] * stress[i][4] * stress[i][4] + stress[i][1] * stress[i][5] * stress[i][5] +
            stress[i][2] * stress[i][3] * stress[i][3]);
            auto node = (*eitr)->nodes[i];
            auto von_mises = std::sqrt(0.5 * ((mStress[i][0] - mStress[i][1]) * (mStress[i][0] - mStress[i][1])
            + (mStress[i][1] - mStress[i][2]) * (mStress[i][1] - mStress[i][2])
            + (mStress[i][2] - mStress[i][0]) * (mStress[i][2] - mStress[i][0])));
            node->ct++;
            node->von_mises += von_mises;
            node->fx += stress[i][0];
            node->fy += stress[i][1];
            node->fz += stress[i][2];
            node->fxy += stress[i][3];
            node->fzy += stress[i][4];
            node->fzx += stress[i][5];
        }
        int i = 0;
    }// for each node in Node Set
    STRAIN_DONE = true;
    STRESS_DONE = true;
    return STRAIN_DONE;
}

 

//========== Solve

// Solve the system, updating the node displacements
// COLLECTIVE

// only call this function if setup done and/or rebuild done

// Solve linear elastic problem
bool vFESolver::Solve(){
    
    idxType ii(0);
    xyzType xx(0), yy(0), zz(0);
    idxType  nodeidx;
    FILE * outFile;
    
    // PC prec;              // preconditioner
    SymetrixSparseMatrix GSM;              // GSM
    Vector sol;              // solution
    Vector rhs;              // rhs (force) vector
    int iters = 0;       // number of iterations
    float norm = 0;       // norm

    printf("In solve \n");
    

    Solver solver;
    
    vFESolver::ComputeGSM(GSM);
    vFESolver::ComputeRHS(rhs);
    
    
    vFESolver::AllocateLocalVec(sol);
    
    
    // KSPSetType(ksp,KSPCG);
    // PCSetType(prec,PCJACOBI);
    // KSPSetTolerances(ksp, TOLERANCE, PETSC_DEFAULT, PETSC_DEFAULT, MAXITER);
    solver.setMaxIteration(MAXITER);
    solver.setTolerance(TOLERANCE);
	// KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    
    
    printf("Starting solver...\n");
    
    auto startT = std::chrono::high_resolution_clock::now();
    
    // KSPSolve(ksp, rhs, sol);
    solver.Solve(GSM,sol,rhs);

    auto endT = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endT - startT);
    
    printf("Done solving!\n");
    printf("SOLVETIME = %.5le \n", static_cast<double>(duration.count()));
    
    // KSPGetSolution(ksp, &sol);
    iters = solver.getIteration();
    norm = solver.getResdualNorm();
    
    // KSPGetIterationNumber(ksp, &iters);
    // KSPGetResidualNorm(ksp, &norm);
    
    
    printf("Converged to %f in %d iterations.\n", norm, iters);
    
    printf("About to print results\n");
    
    Scalar tmpx, tmpy, tmpz;
    int nlocal;
    // VecGetLocalSize(sol, &nlocal);
    
    
    int vstart = 0, vend = 0;
    // VecGetOwnershipRange(sol, &vstart, &vend);
    
    // VecScatter vsctx;
    // VecScatterCreateToZero(sol, &vsctx, &vecout);
    // VecScatterBegin(vsctx, sol, vecout, INSERT_VALUES, SCATTER_FORWARD);
    // VecScatterEnd(vsctx, sol, vecout, INSERT_VALUES, SCATTER_FORWARD);
     
    printf("vstart = %d, vend = %d \n", vstart, vend);

    // VecGetArray(vecout, &solution);
    vecout.reset(sol.size());
    solution = sol.generateScalar();
    
    int vsize;
    vsize = vecout.size();
    printf("vecout size = %d \n", vsize);
    
    // update nodes with final displacements
    NodeSet_it nitr = NodeS.begin();
    Node<xyzType>* nodeptr;

    for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr) {
        nodeptr = (*nitr);
        
        xx = nodeptr->x;
        yy = nodeptr->y;
        zz = nodeptr->z;
        ii = nodeptr->idx;
        // update node displacements
        nodeptr->dx = solution[ii * DOF_3D + 0];
        nodeptr->dy = solution[ii * DOF_3D + 1];
        nodeptr->dz = solution[ii * DOF_3D + 2];
        //printf("%llu %lf %lf %lf \n",ii, solution[ii * DOF_3D + 0], solution[ii * DOF_3D + 1], solution[ii * DOF_3D + 2]);
    }// for each node in Node Set
    
    // VecRestoreArray(vecout, &solution);
    // solution = vecout.generateScalar();
    
    printf("Leaving Solve\n");
    
    SOLVE_DONE = true;
    
    return true;
}// Solve()


//=========== Finish up

bool vFESolver::cleanup()
{
    
    // VecDestroy(&vecout);
    vecout.destroy();
    
    // check if set
    // KSPDestroy(&ksp);
    printf("All done, bye, bye\n");
    
    return true;
}// cleanup()


// Renumber nodes in NodeSet
bool vFESolver::RenumberNodes()
{
    bool OK(true);
    
    try
    {
        idxType ncount(0);
        idxType tmpidx(0);
        NodeSet_it nitr;

        for (nitr = NodeS.begin(); nitr != NodeS.end(); ++nitr)
        {
            tmpidx = (*nitr)->idx;
            (*nitr)->idx = ncount;
            ++ncount;
        }// for each node
    }
    catch (std::exception& e)
    {
        printf("ERROR: could not renumber node\n");
        OK = false;
    }
    return OK;
}

