//
//  VFEModel.cpp
//
//
//

#include "VFEModel.hpp"


//========== CONSTRUCTOR / DESTRUCTOR

VFEModel::VFEModel()
{
}

VFEModel::~VFEModel()
{}

// script command
void VFEModel::create_command_map()
{
    
    // **** Must match enums in Commands.h ****
    // MUST have all options set before LOAD CONSTRAINTS called
    
    // dud command for sake of differentiating between different script formats
    CommandM["SET_SCRIPT_VERSION"] = CMD_SET_SCRIPT_VERSION;
    
    CommandM["SET_VOXEL_SIZE"] = CMD_SET_VOXEL_SIZE;
    CommandM["VOXEL_SIZE"] = CMD_SET_VOXEL_SIZE;
    
    CommandM["LOAD_MATERIALS"] = CMD_LOAD_MATERIALS;
    CommandM["LOAD_MATERIALS_FILE"] = CMD_LOAD_MATERIALS;
    
    CommandM["LOAD_MODEL"] = CMD_LOAD_MODEL;
    CommandM["LOAD_MCTSCAN"] = CMD_LOAD_MODEL;
    
    CommandM["SET_TOLERANCE"] = CMD_SET_TOLERANCE;
    CommandM["TOLERANCE"] = CMD_SET_TOLERANCE;
    
    CommandM["SET_MAX_ITER"] = CMD_SET_MAX_ITER;
    CommandM["MAX_ITER"] = CMD_SET_MAX_ITER;
    
    CommandM["SET_ALGORITHM_FEA"] = CMD_SET_ALGORITHM_FEA;
    CommandM["ALGORITHM_FEA"] = CMD_SET_ALGORITHM_FEA;
  
    // last command before solve is prompted
    CommandM["SELECTION_OF_NODES"] = CMD_SELECTION_OF_NODES; // DOES NOTHING !!
    CommandM["LOAD_CONSTRAINTS"] = CMD_LOAD_CONSTRAINTS;
    CommandM["SELECT_NODE_FILE"] = CMD_LOAD_CONSTRAINTS;
    
    CommandM["SELECT_NODE_3D"] = CMD_SELECT_NODE_3D;
    CommandM["PRESERVE_NODE_3D"] = CMD_PRESERVE_NODE_3D;
    
    CommandM["COMPUTE_SED"] = CMD_COMPUTE_SED; // DOES NOTHING !!
    
    CommandM["SOLVE"] = CMD_SOLVE; // not specified in script
    
    CommandM["PRINT_DISPLACEMENTS"] = CMD_PRINT_DISPLACEMENTS;
    CommandM["PRINT_X"] = CMD_PRINT_DISPLACEMENTS;
    
    CommandM["FINISH"] = CMD_FINISH; // not specified in script
}// create_command_map()


//
// break input string into command and data
void VFEModel::parse_input_line(std::string line)
{
    int i(0);
    int index(0);
    
    size_t cstart = line.find_first_not_of(" \t\n\r");

    
    // only add line if non-white space characters found
    if( cstart != std::string::npos )
    {
        size_t cend = line.find_first_of(" \t\n\r", cstart);
        std::string commandstr = line.substr(cstart, cend-cstart);
       
        // search for command in command map
        if ((CommandM.find(commandstr) != CommandM.end()))
            index = CommandM[ commandstr ];
            else
        {
            printf("ERROR: command %s does not exist, please check script file\n", commandstr.c_str());
            return;
        }
        
        if(cend < line.length())
            InputM[ index ] = line.substr(cend, line.length());
        else
            InputM[ index ] = "";
        
        printf("COMMAND %d:  %s,  DATA: %s \n", index, commandstr.c_str(), (InputM[index]).c_str() );
    }// if non-empty line
    
}// parse_input_line()


// MASTER ONLY
void VFEModel::output_end_flag(bool allOK)
{
    char filename[] = "FECompleted.txt";
    
    FILE * finalFile = fopen(filename,"w");
    
    if(allOK)
        fprintf(finalFile, "solver_done\n");
    else
        fprintf(finalFile, "solver_failed\n");
    
    fclose(finalFile);
}

// MASTER ONLY
void VFEModel::output_input_flag()
{
    char filename[] = "InputCompleted.txt";
    
    FILE * doneFile = fopen(filename,"w");
    
    fprintf(doneFile, "solver_files_input\n");
    
    fclose(doneFile);
}

//========== execute commands

void VFEModel::read_script_file(const char *filename)
{
    startT = MPI_Wtime(); 
    
    bool          OK(true);
    FILE *        scriptFile;
    
    // set default script/model file format
    SCRIPTVERSION = 2;
    
    
    // SET SCRIPTFILE
    strcpy(solver.SCRIPTFILE, filename);
    
    int command(0);
    int totalscript(0);
    InputM_it itr;
    
    // create command map
    create_command_map();
    
    // read commands from script into Input map
    
    // open script file
    std::ifstream scriptF(filename);
    if(scriptF)
    {
        try
        {
            printf("about to read from the script file : %s \n", filename);
            int i(0);
            
            // for each line in script
            for (std::string line; getline(scriptF, line); )
                parse_input_line(line);

            // add finish command
            InputM[CMD_FINISH] = " ";
            
            totalscript = InputM.size();
            itr = InputM.begin();
        }// try this
        catch( std::exception &e )
        {
            printf("ERROR: could not read script file \n");
            std::cerr << "exception caught: " << e.what() << std::endl;
            OK = false;
            
            output_end_flag(OK);
            
            return;
        }// catch this
        scriptF.close();
    }// if file exists
    else
    {
        printf("ERROR: could not find script file \n");
        return;
    }
    
    itr = InputM.begin();
        
    
    
    printf("script read, totalinput = %d \n",totalscript);
    
    int in(0);
    int count(0);
    
    
    // execute commands
    while( OK && (command != CMD_FINISH) )
    {// while OK and not finished
        
        // check if setup done, else read from script input
        {// get command
            
            if( solver.IsSetupDone() && !solver.SOLVE_DONE )
            {
                printf("set up is done, about to solve \n");
                command = CMD_SOLVE;
                // no data string for SOLVE;
                
                //output input file flag
                output_input_flag();
                while (((itr->first!=command)))
                {
                    itr++;
                }
                itr++;
            }
            else
            {
                command = (itr->first);
                strcpy( data, (itr->second).c_str());
                
                ++itr;
            }
        }// if MASTER
        
        
        OK = execute_command(command);
        printf("command (%d) executed and OK = %d\n", command, OK);
        
    }// for each input line
    
    if( !OK ) // if something went wrong with setup
    {
        printf("ERROR: could not continue reading script, something went wrong \n");
        return;
    }
    else
    {
        printf("all OK \n");
    }
    
    endT = MPI_Wtime();
    
    output_end_flag(OK);
    printf("total time = %le \n", endT - startT);
    
}// read_script_file();


// execute commands here
bool VFEModel::execute_command(int command)
{
    //printf("%d: about to execute with command = %d \n", MPIrank, command);
    bool OK(false);
    switch (command)
    {
        
        case CMD_SET_SCRIPT_VERSION : // SET SCRIPT VERSION --
        {
            OK = ( sscanf( data, " %d \n", &SCRIPTVERSION ) == 1 );
            // if not 2, assume default, i.e. previous, formats (default = 1)
            // NOTE: any int that's not 2 indicates default
            if( OK )
            {
                solver.SCRIPTVERSION = SCRIPTVERSION ;
            }
            break;
        }
        case CMD_SET_VOXEL_SIZE : // SET VOXEL SIZE, SCALE FACTOR
        {
            double vsize[4]; // {a, b, c, sf};
            OK = ( sscanf( data, " %lf %lf %lf %lf \n", &vsize[0], &vsize[1], &vsize[2], &vsize[3]) == 4 );
            if (OK)
            {
                OK = solver.SetVoxelSize( vsize[0], vsize[1], vsize[2], vsize[3] );
            }
            break;
        }
        case CMD_LOAD_MATERIALS : // LOAD MATERIALS FROM FILE
        {
            char materialsfn[MAX_FILENAME_LENGTH];
            int materialsfnlen(0);
            OK = ( sscanf( data, " %s \n", &materialsfn ) == 1 );
            materialsfnlen = strlen(materialsfn) + 1;
            if (OK)
            {
                OK = solver.LoadMaterials( materialsfn );
            }
            break;
        }
        case CMD_LOAD_MODEL: // SET VOXEL DIMENSIONS, LOAD MODEL
        {
            xyzType vdims[3];
            xyzType totelems;
            char modelfn[MAX_FILENAME_LENGTH];
            int modelfnlen(0);
            if( SCRIPTVERSION != 2 ) // default
                OK = (sscanf(data, " %d %d %d %s \n", &vdims[0], &vdims[1], &vdims[2], &modelfn) == 4);
            else
                OK = (sscanf(data, " %d %d %d %d %s \n", &vdims[0], &vdims[1], &vdims[2], &totelems, &modelfn) == 5);
                
            modelfnlen = strlen(modelfn) + 1;
            if (OK)
            {
                OK = solver.SetVoxelDimensions( vdims[0], vdims[1], vdims[2]);
                if (OK)
                {
                    if (SCRIPTVERSION == 2)
                    {
                        OK = solver.SetTotElems( totelems );
                    }
                }
            }
            if (OK)
            {
                OK = solver.LoadModel(modelfn);
            }
            break;
        }

        case CMD_SET_TOLERANCE: // TOLERANCE
        {
            double thetolerance;
            OK = (sscanf(data, " %lf \n", &thetolerance) == 1);
            if (OK)
            {
                OK = solver.SetTolerance(thetolerance);
            }
            break;
        }
    	
        case CMD_SET_MAX_ITER : // MAX ITERATION
        {
            long int themaxiter;
            OK = (sscanf(data, " %d \n", &themaxiter) == 1);
            if (OK)
            {
                OK = solver.SetMaxIter( themaxiter );
            }
            break;
        }

        case CMD_SET_ALGORITHM_FEA: // SET SOLVERTYPE & PCTYPE
        {
            char thesolvertype[MAX_FILENAME_LENGTH];
            char thepctype[MAX_FILENAME_LENGTH];
            int solvertypelen;
            int pctypelen;
            int nargs;
            nargs = sscanf(data, " %s  %s \n", &thesolvertype, &thepctype);

            if (nargs == 1)
            {
                strcpy(thepctype, "");
                OK = true;
            }
            if (nargs == 2)
                OK = true;

            solvertypelen = strlen(thesolvertype);
            pctypelen = strlen(thepctype);

            if (OK)
            {
                OK = solver.SetAlgorithmFEA(thesolvertype, thepctype);
            }
            break;
        }

        case CMD_SELECTION_OF_NODES: // SELECTION OF NODES 
        {
            OK = true;
            break;
        }

        case CMD_LOAD_CONSTRAINTS: // LOAD CONSTRAINTS (FORCES, PRESERVED NODES)
        {
            char constraintsfn[MAX_FILENAME_LENGTH];
            int constraintsfnlen(0);
            OK = (sscanf(data, " %s \n", &constraintsfn) == 1);
            constraintsfnlen = strlen(constraintsfn) + 1;
            if (OK)
            {
                OK = solver.LoadConstraints(constraintsfn);
            }
            break;
        }
        case CMD_TRANSFER_TO_NASTRAN:
        {
            char NastranFN[MAX_FILENAME_LENGTH];
            int Nastranfnlen(0);
            OK = (sscanf(data, " %s \n", &NastranFN) == 1);
            Nastranfnlen = strlen(NastranFN) + 1;
            if (OK)
            {
                OK = solver.toNASTRAN(NastranFN);
            }
            break;
        }
        case  CMD_SELECT_NODE_3D: //
        {
            OK = true;
            break;
        }

        case CMD_PRESERVE_NODE_3D: //
        {
            OK = true;
            break;
        }
    	
        case CMD_COMPUTE_SED : //
        {
            OK = true;
            break;
        }
        
        case CMD_SOLVE : // SOLVE
        {
            OK = (!solver.Solve()); // PetscErrorCode
            break;
        }
            
        case CMD_PRINT_DISPLACEMENTS : // PRINT DISPLACEMENTS
        {
            // only MASTER process prints displacements (for now)
            if (solver.SOLVE_DONE)
            {
                char outputfn[MAX_FILENAME_LENGTH];
                OK = (sscanf(data, " %s \n", &outputfn) == 1);
                if (OK)
                    OK = solver.PrintDisplacements(outputfn);
            
            }// if solve done
            else
            {
                printf( "%d: ERROR: cannot print displacements before solving! \n");
                OK = false;
            }
            
            break;
        }

        case CMD_FINISH : // FINISHED READING SCRIPT
        {
            printf( "Finished reading script \n");
            // cleanup here
            solver.FINISH_DONE = true;
            OK =  true;
            
            break;
        }
        default:
            printf("ERROR: command %d not recognised\n", command);
    }
    
    return OK;
}// execute_command()