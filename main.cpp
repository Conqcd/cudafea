#include "VFEModel.hpp"

using namespace std;

int main(int argc, char* argv[])
{   
    // Start up

    VFEModel *vox;
    vox = new VFEModel();
    
    printf("calling read_script_file \n");
    vox -> read_script_file(argv[1]);
    
    //clean up and shut down
    delete vox;
    vox = 0;
    
    return 0;
}
