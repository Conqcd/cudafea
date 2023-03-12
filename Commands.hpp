//
//  Commands.h
//

#ifndef COMMANDS_H
#define COMMANDS_H


#define MAX_COMMAND_LENGTH 256
#define MAX_DATA_LENGTH 256
#define MAX_FILENAME_LENGTH 256

// Rank of the master process
enum { MASTER = 0 };

//======================== COMMANDS =====================================

enum {  CMD_SET_SCRIPT_VERSION, 
        CMD_SET_VOXEL_SIZE, 
        CMD_LOAD_MATERIALS, 
        CMD_LOAD_MODEL, 
        CMD_SET_TOLERANCE, 
        CMD_SET_MAX_ITER, 
        CMD_SET_ALGORITHM_FEA, 
        CMD_SELECTION_OF_NODES, 
        CMD_GENERATE_CONSTRAINTS, 
        CMD_LOAD_CONSTRAINTS, 
        CMD_TRANSFER_TO_NASTRAN,
        CMD_SELECT_NODE_3D, 
        CMD_PRESERVE_NODE_3D, 
        CMD_COMPUTE_SED, 
        CMD_SOLVE,
		CMD_GET_FORCE,
        CMD_PRINT_DISPLACEMENTS, 
        CMD_PRINT_STRAIN, 
        CMD_PRINT_STRESS, 
        CMD_FINISH };


#endif /* COMMANDS_H */
