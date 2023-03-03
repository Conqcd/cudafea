//
//  VFEMDOEL.h
//

#ifndef VFEMODEL_H
#define VFEMODEL_H

#include "vFESolver.hpp"

class VFEModel{

public:
    
    int         SCRIPTVERSION;
    
    char data[MAX_DATA_LENGTH];
    
    vFESolver   solver;
    
    typedef        std::map<std::string, int> MapStringInt;  // map of commands
    typedef        std::map<int, std::string> MapIntString;  // map of commands
    
    MapStringInt   CommandM; // mapping commands to long ints
    MapIntString   InputM; // mapping commands to long ints
    
    typedef        MapStringInt::iterator           CommandM_it;
    typedef        MapIntString::iterator             InputM_it;
    
    // measure time
    double startT, endT;
    
    VFEModel();
    ~VFEModel();

    // create map < command string , command enum>
    void create_command_map();
    
    // MASTER ONLY
    void read_script_file(const char * filename);
    
    // parse line from script --> command, data
    void parse_input_line(std::string line);
    
    // COLLECTIVE
    bool execute_command(int command);
    
    
    // output final flag file to show solver has completed 
    void output_end_flag(bool allOK);
    
   
    // output flag file once input has been read in
    void output_input_flag();
};



#endif /* VFEMODEL_H */
