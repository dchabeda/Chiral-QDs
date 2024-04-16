#include "fd.h"
#define MAX_LINE_LENGTH 200

int main(int argc, char *argv[]){
    int chkIdx = atoi(argv[1]);
    printf("The chkIdx = %ld\n", chkIdx);

    if (chkIdx == 1){
        // Restart the job after obtaining all the filtered states
        return 0;
    }

    if (chkIdx == 2){
        // Restart job from after the full Hmat diagonalization. psi.dat has been printed.
        FILE *pf;
        char line[MAX_LINE_LENGTH]; //Arbitrary max line length
        int ieof;

        //Make sure that the psi.dat finished writing out to disk
        pf = fopen("run.dat", "r");

        if (pf == NULL){
            printf("No run.dat in directory! Checkpoint restart failed.");
            return 1;
        }
        while (fgets(line, MAX_LINE_LENGTH, pf)){
            printf(line);
        }

        return 0;
    }
    printf("Invalid chkIdx. Checkpoint restart failed.\n");
    return 1;
}