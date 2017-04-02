#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>
#include <assert.h>
#include <chrono>
#include <ctime>
#include <vector>
#include <stddef.h>
#include <sys/time.h>
#include <unistd.h>


#define G3D(V,X,Y,Z)  V[(Z) * ((TAM_X) * (TAM_Y)) + (Y) * (TAM_X) + (X)]
#define S3D(V,X,Y,Z,S)  V[(Z) * ((TAM_X) * (TAM_Y)) + (Y) * TAM_X + (X)]=S
#define CREATEM3D(DIMX,DIMY,DIMZ) std::vector<double>((DIMX)*(DIMY)*(DIMZ))
#define VECTOR3D std::vector<double>

int TAM_X = 90;
int TAM_Z = 90;
int TAM_Y = 122;

double CM = 1E+8;


using namespace std;



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
    
    osmc << "ORIGIN    45.000   50.000   45.000" << endl;
    osmc << "SPACING    1.000   1.000   1.000" << endl;
    
    osmc << "POINT_DATA   " << numDataPoints << endl;
    osmc << "SCALARS scalars double" << endl;
    osmc << "LOOKUP_TABLE default" << endl;
    
    for (int i = 0; i < mat.size();i++)
        osmc << mat[i] << endl;
    
    
    fbmc.close();
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

int main()
{
    


    
    VECTOR3D cerebro = CREATEM3D(TAM_X,TAM_Y,TAM_Z);
    VECTOR3D difusion = CREATEM3D(TAM_X,TAM_Y,TAM_Z);
    
    ReadDifussionData("./Cerebro.csv", TAM_X, TAM_Y, TAM_Z, 45, 50, 45, difusion);
    
    
    dumpMatrixToVtk(difusion, "difusion");
    
    printf ("Preprocessing difusion Matrix\n");
    TransformDifusion(difusion);
    
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
    dumpMatrixToVtk(cerebro, "cerebro_0");
   

    
    double dx = 1;double dy=1;double dz=1;
    double dx2 = dx * dx;double dy2 = dy * dy;double dz2 = dz * dz;
    double dt = 0.003;
    
    VECTOR3D DDX = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
    VECTOR3D DDY = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
    VECTOR3D DDZ = CREATEM3D(TAM_X, TAM_Y, TAM_Z);
    
    GetDifussionDerivateX(difusion, DDX, dx);
    GetDifussionDerivateY(difusion, DDY, dy);
    GetDifussionDerivateZ(difusion, DDZ, dz);
    
    dumpMatrixToVtk(DDX, "DDX");
    
    VECTOR3D& D = difusion;
    
    double epsilon = 0.01;
    
    printf("Concentración total: %f\n", getConcentracionTotal(cerebro));
    
    double t1,t2,elapsed;
    struct timeval tp;
    int rtn;

    rtn=gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;    
    
    for (int t = 1; t < 100; t++){
        
        //cout << GetLocalTime() << endl;
        
        //printf ("Starting iteration %u....\n", t);
        /*
        for (int k = 1; k < TAM_Z-2; k++)
        {
                for (int j = 1; j < TAM_Y-2; j++)
                {
                    for (int i = 1; i < TAM_X-2; i++)
*/
        for (int i = 1; i < TAM_X-2; i++)
        {
                for (int j = 1; j < TAM_Y-2; j++)
                {
                    for (int k = 1; k < TAM_Z-2; k++)

                    {
                         double c = G3D(cerebro, i, j, k);
                         double d = G3D(D, i, j, k);
                         double c_plus_i =G3D(cerebro,i+1,j,k);
                         double c_minus_i =G3D(cerebro,i-1,j,k);
                         double c_plus_j =G3D(cerebro,i,j+1,k);
                         double c_minus_j =G3D(cerebro,i,j-1,k);
                         double c_plus_k =G3D(cerebro,i,j,k+1);
                         double c_minus_k =G3D(cerebro,i,j,k-1);
                        
                        
                         double c_new =
                        
                        c
                        
                        +
                        (
                         G3D(DDX, i, j, k) * (c_minus_i - c_plus_i) / (2 * dx) + d * (c_plus_i - 2*c + c_minus_i) / dx2
                         +
                         G3D(DDY, i, j, k) * (c_minus_j - c_plus_j) / (2 * dy) + d * (c_plus_j - 2*c + c_minus_j) / dy2
                         +
                         G3D(DDZ, i, j, k) * (c_minus_k - c_plus_k) / (2 * dz) + d * (c_plus_k - 2*c + c_minus_k) / dz2
                         +
                         RO(d) * c * log(CM / sqrt(c*c + epsilon))
                         
                         - R(c)
                         )
                        
                        * dt;
 
                        if (c_new > 1E8)
                            c_new = 1E8;
                        if (c_new < 0) c_new = 0;
                        if (c_new < 1e-300) c_new = 0;

                        
                        if (c_new != c_new)
                            printf ("T: %u", t);
                
                        
                        S3D(cerebro,i,j,k,c_new);
                        
                        
                    }
                }
        }
    

        
        
        if (t % 500 == 0){
            dumpMatrixToVtk(cerebro, "cerebro_" + std::to_string(t));
        }
        
        
        /*
        // Neuman
        for (int i = 1; i < TAM_X-2; i++)
        {
            for (int j = 1; j < TAM_Y-2; j++){
                    S3D(cerebro,i,j,0,G3D(cerebro,i,j,1));
                    S3D(cerebro,i,j,TAM_Z-1,G3D(cerebro,i,j,TAM_Z-2));

            }
        }
        
        for (int i = 1; i < TAM_X-2; i++)
        {
            for (int k = 1; k < TAM_Z-2; k++){
                S3D(cerebro,i,0,k,G3D(cerebro,i,1,k));
                S3D(cerebro,i,TAM_Y-1,k,G3D(cerebro,i,TAM_Y-1,k));
                
            }
        }

        for (int j = 1; j < TAM_Y-2; j++)
        {
            for (int k = 1; k < TAM_Z-2; k++){
                S3D(cerebro,0,j,k,G3D(cerebro,1,j,k));
                S3D(cerebro,TAM_X-1,j,k,G3D(cerebro,TAM_X-2,j,k));
                
            }
        }*/
        
        
    }
 
         rtn=gettimeofday(&tp, NULL);
        t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
        elapsed=t2-t1;

        printf("Tiempo empleado: %g\n",elapsed);   
    
    printf("Concentración total: %f\n", getConcentracionTotal(cerebro));

    dumpMatrixToVtk(cerebro, "cerebro_out");

    
}

