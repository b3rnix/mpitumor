#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>
#include <assert.h>
#include <chrono>
#include <ctime>
#include <vector>
#include <mpi.h>
#include <stdarg.h>

#include <stddef.h>
#include <sys/time.h>
#include <unistd.h>


#define G3D(V,X,Y,Z)  V[(Z) * ((TAM_X) * (TAM_Y)) + (Y) * (TAM_X) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((TAM_X) * (TAM_Y)) + (Y) * TAM_X + (X)]=S
#define CREATEM3D(DIMX,DIMY,DIMZ) std::vector<double>((DIMX)*(DIMY)*(DIMZ))
#define VECTOR3D std::vector<double>

int TAM_X = 90;
int TAM_Y = 90;
int TAM_Z = 122;

double CM = 1E+8;


using namespace std;



/*
std::string GetLocalTime() {
    auto now(std::chrono::system_clock::now());
    auto seconds_since_epoch(
                             std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()));
    
    // Construct time_t using 'seconds_since_epoch' rather than 'now' since it is
    // implementation-defined whether the value is rounded or truncated.
    std::time_t now_t(
                      std::chrono::system_clock::to_time_t(
                                                           std::chrono::system_clock::time_point(seconds_since_epoch)));
    
    char temp[10];
    if (!std::strftime(temp, 10, "%H:%M:%S.", std::localtime(&now_t)))
        return "";
    
    return std::string(temp) +
    std::to_string((now.time_since_epoch() - seconds_since_epoch).count());
}
*/





// En realidad tendria que chequear si es materia gris/blanca o no
double RO(double value){
    if ( value >0)
        return 0.107;
    else
        return 0;
}

double R(double c){
    double epsilon = 0.01;
    double alpha = 0.036;
    double beta = 0.0036;
    double dose = 1.8;
    double s = exp(-(alpha * dose) -(beta * dose*dose));
    return (1-s) * c * (double)log(CM / sqrt(c*c + epsilon));
}


void GetDifussionDerivateX(VECTOR3D &D, VECTOR3D &DD, double delta){
    int dimx =TAM_X;
    int dimy =TAM_Y;
    int dimz =TAM_Z;
    
    
    for (int k = 1; k < dimz - 1;k++){
        for (int j = 1; j < dimy - 1;j++){
            for (int i = 1; i < dimx - 1;i++){
                S3D(DD, i, j, k, (G3D(D,i + 1, j, k) - G3D(D,i - 1, j, k)) / (2 * delta));
            }
        }
    }

}

void GetDifussionDerivateY(VECTOR3D &D, VECTOR3D &DD,double delta){
    int dimx =TAM_X;
    int dimy =TAM_Y;
    int dimz =TAM_Z;
    
    
    for (int k = 1; k < dimz - 1;k++){
        for (int j = 1; j < dimy - 1;j++){
            for (int i = 1; i < dimx - 1;i++){
                S3D(DD, i, j, k, (G3D(D,i, j+1, k) - G3D(D,i, j-1, k)) / (2 * delta));
            }
        }
    }
}

void GetDifussionDerivateZ(VECTOR3D &D, VECTOR3D &DD , double delta){
    int dimx =TAM_X;
    int dimy =TAM_Y;
    int dimz =TAM_Z;
    
    
    for (int k = 1; k < dimz - 1;k++){
        for (int j = 1; j < dimy - 1;j++){
            for (int i = 1; i < dimx - 1;i++){
                S3D(DD, i, j, k, (G3D(D,i, j, k+1) - G3D(D,i, j, k-1)) / (2 * delta));
            }
        }
    }
 
}

void dumpMatrixToVtk(VECTOR3D &mat, string fileId){
    
    /*
    cout << "Dumping to VTK..." << endl;
    
    long numDataPoints = TAM_X * TAM_Y * TAM_Z;
    
    std::filebuf fbmc;
    fbmc.open ("./" + fileId + ".vtk",std::ios::out);
    std::ostream osmc(&fbmc);
    
    osmc << "# vtk DataFile Version 2.0" << endl;
    osmc << "Comment goes here" << endl;
    osmc << "ASCII" << endl;
    
    osmc << "DATASET STRUCTURED_POINTS" << endl;
    osmc << "DIMENSIONS    " <<  TAM_X << " " << TAM_Y << " " << TAM_Z << endl;
    
    osmc << "ORIGIN    0.000   0.000   0.000" << endl;
    osmc << "SPACING    1.000   1.000   1.000" << endl;
    
    osmc << "POINT_DATA   " << numDataPoints << endl;
    osmc << "SCALARS scalars double" << endl;
    osmc << "LOOKUP_TABLE default" << endl;
    
    for (int i = 0; i < mat.size();i++)
        osmc << mat[i] << endl;
    
    
    fbmc.close();*/
}

void blurMatrix3d(VECTOR3D& m, int i, int j, int k){
    
    for (int x = 1; x < i; x++){
        for (int y = 1; y < j; y++){
            for (int z = 1; z < k; z++){
                double val = (G3D(m,x,y,z)+G3D(m,x+1,y,z)+G3D(m,x-1,y,z)+G3D(m,x,y+1,z)+G3D(m,x,y-1,z)+G3D(m,x,y,z+1)+G3D(m,x,y,z-1))/8.0;
                S3D( m, x,  y,  z,  val);
            }
        }
    }
    
}

void TransformDifusion(VECTOR3D &src){
    double transformed;
    int tamx = TAM_X;
    int tamy = TAM_Y;
    int tamz = TAM_Z;
    
    for (int k = 1; k < tamz;k++){
        for (int j = 1; j < tamy;j++){
            for (int i = 1; i < tamx;i++){
                double value = G3D(src,i, j, k);
                if ( value >= 110 && value < 255)
                    transformed = 0.255;
                else if (value <= 110 && value >= 55)
                    transformed = 0.051;
                else
                    transformed = 0;
            
                S3D(src, i, j, k, transformed);
            }
        }
    }
    
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);
    blurMatrix3d(src,tamx -1, tamy -1, tamz - 1);

    
}



void ReadDifussionData(string dataFile, int tamX, int tamY, int tamZ, int originX, int originY, int originZ, VECTOR3D &difusionMat){

    printf ("Reading difussion data file %s, section (%u,%u,%u),(%u,%u,%u)\n", dataFile.c_str(), originX, originY, originZ, originX + tamX, originY+tamY, originZ+tamZ);

    
    ifstream file ( dataFile );
    string value;
    getline ( file, value);
    int x = 0;
    int y = 0;
    int z = 0;
    
 
    while ( x < originX + tamX){
                getline ( file, value, ',' );
                x = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );
        
                y = atoi(value.c_str()) - 1;
                getline ( file, value, ',' );
        
                z = atoi(value.c_str()) - 1 ;
                getline ( file, value );
        
                double s = atof(value.c_str());
        
                //printf("%u,%u,%u\n", x,y,z);
                if (originX <= x && x - originX < tamX  && originY <= y && y - originY < tamY && originZ <= z && z - originZ < tamZ)
                    S3D(difusionMat,x - originX,y - originY,z - originZ,s);
            
    }
    printf ("Difussion data read OK\n");
    //S3D(difusionMat,0,0,0,66);

}

double getConcentracionTotal(VECTOR3D mat){
    double ret = 0;
    for (int i = 1; i < TAM_X - 1; i++)
    {
        for (int j = 1; j < TAM_Y - 1; j++)
        {
            for (int k = 1; k < TAM_Z -1; k++){
                ret+= G3D(mat, i, j, k);
            }
        }
    }
    return ret;
}


void screenMessage(const char *fmt, ...) {
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_rank == 0){
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
}

void doHaloTransfer(VECTOR3D &slice, int rank, int numProcs, int workingSliceSize){
     //printf("[%u] Entering Halo\n", rank);
    
     if (rank % 2 == 0){
        //Z0-Z1-Z2-Z3-Z4  Z3-Z4-Z5-Z6-Z7  Z6-Z7-Z8-Z9-ZA  Z9-ZA-ZB-ZC ....
        // SEND Z4(0) <- Z4(1)
        if (rank > 0){
            //printf("[%u] [Send-Left] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Send(&slice[TAM_X  * TAM_Y], TAM_X  * TAM_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            //printf("[%u] [Send-Left] END Halo size: %u\n", rank, workingSliceSize);

        }

        // SEND Z3(1) <- Z3(0)

        if (rank < numProcs - 1){
            //printf("[%u] [Send-Right] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Send(&slice[(TAM_X  * TAM_Y) * (workingSliceSize - 2)], TAM_X  * TAM_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            //printf("[%u] [Send-Right] END Halo size: %u\n", rank, workingSliceSize);

        }

        MPI_Barrier(MPI_COMM_WORLD);


        // RECV Z3(1) <- Z3(0)
        if (rank > 0){

            //printf("[%u] [Recv-Left] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Recv(&slice[0], TAM_X  * TAM_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("[%u] [Recv-Left] END Halo size: %u\n", rank, workingSliceSize);

        }

        // RECV Z4(0) <- Z4(1)
        if (rank < numProcs - 1){
            // printf("[%u] [Recv-Right] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Recv(&slice[(TAM_X  * TAM_Y) * (workingSliceSize-1)], TAM_X  * TAM_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("[%u] [Recv-Right] END Halo size: %u\n", rank, workingSliceSize);

        }
    }else{
            //Z0-Z1-Z2-Z3-Z4  Z3-Z4-Z5-Z6-Z7  Z6-Z7-Z8-Z9-ZA  Z9-ZA-ZB-ZC ....
        // SEND Z4(0) <- Z4(1)
        if (rank > 0){
            //printf("[%u] [Send-Left] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Send(&slice[TAM_X  * TAM_Y], TAM_X  * TAM_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            //printf("[%u] [Send-Left] END Halo size: %u\n", rank, workingSliceSize);

        }

        // SEND Z3(1) <- Z3(0)

        if (rank < numProcs - 1){
            //printf("[%u] [Send-Right] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Send(&slice[(TAM_X  * TAM_Y) * (workingSliceSize - 2)], TAM_X  * TAM_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            //printf("[%u] [Send-Right] END Halo size: %u\n", rank, workingSliceSize);

        }

        MPI_Barrier(MPI_COMM_WORLD);


        // RECV Z3(1) <- Z3(0)
        if (rank > 0){

            //printf("[%u] [Recv-Left] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Recv(&slice[0], TAM_X  * TAM_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("[%u] [Recv-Left] END Halo size: %u\n", rank, workingSliceSize);

        }

        // RECV Z4(0) <- Z4(1)
        if (rank < numProcs - 1){
            // printf("[%u] [Recv-Right] BEGIN Halo size: %u\n", rank, workingSliceSize);
            MPI_Recv(&slice[(TAM_X  * TAM_Y) * (workingSliceSize-1)], TAM_X  * TAM_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("[%u] [Recv-Right] END Halo size: %u\n", rank, workingSliceSize);

        }

    }
     //printf("[%u] Exit Halo\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);

}

int main()
{
    
      // Initialize the MPI environment
      MPI_Init(NULL, NULL);
      // Find out rank, size
      int world_rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      int world_size;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    double dx = 1;double dy=1;double dz=1;
    double dx2 = dx * dx;double dy2 = dy * dy;double dz2 = dz * dz;
    double dt = 0.003;

        VECTOR3D DDX ;
        VECTOR3D DDY ;
        VECTOR3D DDZ ;    
        VECTOR3D cerebro ; 
        VECTOR3D D ;
    if (world_rank == 0) {
        cerebro  = CREATEM3D(TAM_X,TAM_Y,TAM_Z);
        D = CREATEM3D(TAM_X,TAM_Y,TAM_Z);
        ReadDifussionData("./Cerebro.csv", TAM_X, TAM_Y, TAM_Z, 45, 50, 45, D);
        dumpMatrixToVtk(D, "difusion");
        printf ("Preprocessing diffusion Matrix\n");
        TransformDifusion(D);
        // Inicializacion Cerebro
        double radioTumor = 15.0;
        double posTumorX = TAM_X / 2;
        double posTumorY = TAM_Y / 2;
        double posTumorZ = TAM_Z / 2;
        printf ("Creating initial Cerebrus....");

        for (int x = 1; x < TAM_X; x++){
            for (int y = 1; y < TAM_Y; y++){
                for (int z = 1; z < TAM_Z; z++){
                    double norm = sqrt(pow(posTumorX - x,2) +  pow(posTumorY - y,2) +  pow(posTumorZ - z,2));
                    if (norm < radioTumor)
                        S3D(cerebro,x,y,z,(1E8 / 2));
                    else
                        S3D(cerebro,x,y,z,  0.0);
                    
                }
            }
        }
        //dumpMatrixToVtk(cerebro, "cerebro0");
       

        DDX = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
        DDY = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
        DDZ = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
        
        GetDifussionDerivateX(D, DDX, dx);
        GetDifussionDerivateY(D, DDY, dy);
        GetDifussionDerivateZ(D, DDZ, dz);
        
        printf("Concentración total: %f\n", getConcentracionTotal(cerebro));
    }

    MPI_Barrier(MPI_COMM_WORLD);

    double epsilon = 0.01;

    
    // Capas a procesar (el contenido )
    int chunkSize = (TAM_Z / world_size) ;
    int scatterSize = chunkSize - 1;
    int workingSize = chunkSize + 1;
    VECTOR3D cerebro_slice  = CREATEM3D(TAM_X, TAM_Y, workingSize);
    VECTOR3D D_slice        = CREATEM3D(TAM_X, TAM_Y, workingSize);
    VECTOR3D DDX_slice      = CREATEM3D(TAM_X, TAM_Y, workingSize);
    VECTOR3D DDY_slice      = CREATEM3D(TAM_X, TAM_Y, workingSize);
    VECTOR3D DDZ_slice      = CREATEM3D(TAM_X, TAM_Y, workingSize);

    screenMessage("Procs: %u, chunkSize: %u, X_SIZE: %u, Y_SIZE: %u, Z_SIZE: %u\n", world_size, chunkSize, TAM_X, TAM_Y, TAM_Z);

    MPI_Barrier(MPI_COMM_WORLD);
        

    // Paso capas internas a los procesos (arranco de Z1)
    screenMessage("Scatering cerebro\n");
    MPI_Scatter(&cerebro[TAM_X*TAM_Y],scatterSize * (TAM_X*TAM_Y),MPI_DOUBLE,&cerebro_slice[TAM_X*TAM_Y],scatterSize * (TAM_X*TAM_Y),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    /*
    screenMessage("Scatering D\n");
    MPI_Scatter(&D[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y),MPI_DOUBLE,&D_slice[TAM_X*TAM_Y],scatterSize * (TAM_X*TAM_Y),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    screenMessage("Scatering DDX\n");
    MPI_Scatter(&DDX[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDX_slice[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    screenMessage("Scatering DDY\n");
    MPI_Scatter(&DDY[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDY_slice[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);

    screenMessage("Scatering DDZ\n");
    MPI_Scatter(&DDZ[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDZ_slice[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    */

    MPI_Scatter(&D[0],chunkSize* (TAM_X*TAM_Y),MPI_DOUBLE,&D_slice[0],chunkSize * (TAM_X*TAM_Y),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    screenMessage("Scatering DDX\n");
    MPI_Scatter(&DDX[0],chunkSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDX_slice[0],chunkSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    screenMessage("Scatering DDY\n");
    MPI_Scatter(&DDY[0],chunkSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDY_slice[0],chunkSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);

    screenMessage("Scatering DDZ\n");
    MPI_Scatter(&DDZ[0],chunkSize* (TAM_X*TAM_Y),MPI_DOUBLE,&DDZ_slice[0],chunkSize* (TAM_X*TAM_Y) ,MPI_DOUBLE,0,MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);

    double t1,t2,elapsed;
    struct timeval tp;
    int rtn;

    if (world_rank == 0){
        rtn=gettimeofday(&tp, NULL);
        t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;    
    }

    for (int t = 1; t < 100; t++){

        //screenMessage("Waiting for processes to border exchange...\n");
        MPI_Barrier(MPI_COMM_WORLD);
        // Intercambio de bordes Cn entre procesos vecinos
        //screenMessage("Starting border exchange...\n");
        doHaloTransfer(cerebro_slice, world_rank, world_size, workingSize);
        MPI_Barrier(MPI_COMM_WORLD);

        
        
        //if (world_rank == 0){
            //cout << GetLocalTime() << endl;
            //printf ("Starting iteration %u....\n", t);
            //printf ("|");
        //}

        for (int k = 1; k < chunkSize; k++)
        {
                for (int j = 1; j < TAM_Y-2; j++)
                {
                    for (int i = 1; i < TAM_X-2; i++)
                    {
                         double c = G3D(cerebro_slice, i, j, k);
                         double d = G3D(D_slice, i, j, k);
                         double c_plus_i =G3D(cerebro_slice,i+1,j,k);
                         double c_minus_i =G3D(cerebro_slice,i-1,j,k);
                         double c_plus_j =G3D(cerebro_slice,i,j+1,k);
                         double c_minus_j =G3D(cerebro_slice,i,j-1,k);
                         double c_plus_k =G3D(cerebro_slice,i,j,k+1);
                         double c_minus_k =G3D(cerebro_slice,i,j,k-1);
                        
                        
                         double c_new =
                        
                        c
                        
                        +
                        (
                         G3D(DDX_slice, i, j, k) * (c_minus_i - c_plus_i) / (2 * dx) + d * (c_plus_i - 2*c + c_minus_i) / dx2
                         +
                         G3D(DDY_slice, i, j, k) * (c_minus_j - c_plus_j) / (2 * dy) + d * (c_plus_j - 2*c + c_minus_j) / dy2
                         +
                         G3D(DDZ_slice, i, j, k) * (c_minus_k - c_plus_k) / (2 * dz) + d * (c_plus_k - 2*c + c_minus_k) / dz2
                         +
                         RO(d) * c * log(CM / sqrt(c*c + epsilon))
                         
                         - R(c)
                         )
                        
                        * dt;
 
                        if (c_new > 1E8)
                            c_new = 1E8;
                        if (c_new < 0) c_new = 0;
                        if (c_new < 1e-300) c_new = 0;
                        
                        S3D(cerebro_slice,i,j,k,c_new);
                        //S3D(cerebro_slice,i,j,k,1000);
                        
                    }
                }
        }
        
    }

    MPI_Barrier(MPI_COMM_WORLD);

    screenMessage("End of loop\n");
    
    MPI_Gather(&cerebro_slice[TAM_X*TAM_Y],scatterSize * (TAM_X*TAM_Y),MPI_DOUBLE,&cerebro[TAM_X*TAM_Y],scatterSize* (TAM_X*TAM_Y),MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    screenMessage("End of data gathering\n");

    
    if (world_rank == 0){
        rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;

        printf("Tiempo empleado: %g\n",elapsed);

        screenMessage("Calculating...\n");
        printf("Concentración total: %f\n", getConcentracionTotal(cerebro));
        dumpMatrixToVtk(cerebro, "cerebro_out");
    }
     MPI_Finalize();    
     return(0);

    
}

