/*******************************************************************************
*
* ants_map
* December 2018
* Takes point sets as input
*
********************************************************************************
*/

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MINX 293706
#define MAXX 564106
#define MINY 6803848
#define MAXY 7076948
#define PI 3.14159265358979323846264338327950288419716939937510582
#define MAXNESTS 50000

#define XOFF 480000
#define YOFF 6940000

#define MAXPOLY 100

int SAMPLESIZE; // after burn in
int MINLEN; //Parameter for chi-shapes
int GRIDSIZE = 100;

char nestCoordsFileName[100];
char outputFileName[100];


/*
* Defines a wrapper class for a 2D array, implemented as a 1D array
*/
template <typename T>
class Array2D{
    // private data members
private:
    int rowNum;
    int columnNum;
    T* Matrix;

public:
    /*
    * Constructor
    */
    Array2D(int row, int column){
        if (row > 0 && column > 0){
            rowNum = row;
            columnNum = column;
            Matrix = new T[rowNum *columnNum + columnNum];
        }

    }

    /*
    * Deconstructor
    */
    ~Array2D(){
        delete[] Matrix;
    }

    /*
    * Returns pointer to item at Matrx[x][y][0]
    * i.e, an "array" of size depthNum
    */
    T* operator () (int x){
        return &Matrix[x*columnNum + 0];
    }

    /*
    * Returns item at Matrix[x][y][z]
    */
    const T& operator () (int x, int y) const{
        return Matrix[x*columnNum + y];
    }
    T& operator () (int x, int y){
        return Matrix[x*columnNum + y];
    }

    int getRowNum(){
        return rowNum;
    }

    int getColumnNum(){
        return columnNum;
    }

};


struct edge{
  int p1; //End point 1
  int p2; //End point 2
  int t1; //occurence of edge in indexlist
};


void handleargs(int argc, char **argv){
    if (argc != 5){
        printf("Invalid commandline.\n");
        exit(1);
    }
    //since nestCoordsFileName[100], the programs terminates if file name length>100
    if (strlen(argv[1]) > 100){
        printf("Nest coords file name >100.\n");
        exit(1);
    }
    //Copy string, from argv[1] into the array nestCordsFileName
    strcpy(nestCoordsFileName, argv[1]);
    
    //Read formatted data from string
    sscanf(argv[2], "%d", &SAMPLESIZE);
    
    sscanf(argv[3], "%d", &MINLEN);
    
    sscanf(argv[4], "%d", &GRIDSIZE);
    
    sprintf(outputFileName,"%s.L%d.polygon",nestCoordsFileName,MINLEN);
}


double *comparelen;
int comparelength(const void *e1,const void *e2)
{
  //Sort by decreasing length, then by increasing p1, then increasing p2.
  struct edge *edge1 = (struct edge*)e1;
  struct edge *edge2 = (struct edge*)e2;
  
  double len1 = comparelen[edge1->t1];
  double len2 = comparelen[edge2->t1];
  
  if(len1 > len2) return -1;
  else if(len1 < len2) return 1;
  else if(edge1->p1 == edge2->p1) return edge1->p2 - edge2->p2;
  else return edge1->p1 - edge2->p1;
}


double *boundarypoints(float *points,int numpoints,int *numboundarypoints,int *numpoly,int minlen)
{
  unsigned int *indexlist;
  int numTriangleVertices;
  int vnum,next;
  int *TriangleCount;
  struct edge *edgelist;
  double *len;
  double x1,y1,x2,y2;
  int *eext;
  int *partner;
  int *removed;
  int edge,edgeA,edgeB;
  int numused;
  double *boundary;


  unsigned int *BuildTriangleIndexList (void *pointList, float factor, int numberOfInputPoints, int numDimensions, int clockwise, int *numTriangleVertices );
  

  indexlist = BuildTriangleIndexList (points, 1.0, numpoints, 2, 1, &numTriangleVertices);
    
  //Count triangles for each vertex
  TriangleCount = (int*)malloc(numpoints*sizeof(int));
  for(vnum=0;vnum<numpoints;vnum++){
    TriangleCount[vnum] = 0;
  }
  for(vnum=0;vnum<numTriangleVertices;vnum++){
    TriangleCount[indexlist[vnum]]++;
  }

  //Build list of edges
  edgelist = (struct edge*)malloc(numTriangleVertices*sizeof(struct edge));
  len = (double*)malloc(numTriangleVertices*sizeof(double));
  for(vnum=0; vnum < numTriangleVertices; vnum++){
    edgelist[vnum].p1 = indexlist[vnum];
    if(vnum%3 == 2) edgelist[vnum].p2 = indexlist[vnum-2];
    else edgelist[vnum].p2 = indexlist[vnum+1];
    if(edgelist[vnum].p2 < edgelist[vnum].p1){
      //Swap vertices
      edgelist[vnum].p1 = edgelist[vnum].p2;
      edgelist[vnum].p2 = indexlist[vnum];
    }
    edgelist[vnum].t1 = vnum;
    
    x1 = points[2*(edgelist[vnum].p1)];
    y1 = points[2*(edgelist[vnum].p1)+1];
    x2 = points[2*(edgelist[vnum].p2)];
    y2 = points[2*(edgelist[vnum].p2)+1];
    len[vnum] = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  }
  
  //sort by length
  comparelen = len;
  qsort(edgelist,numTriangleVertices,sizeof(struct edge),&comparelength);
  
  //identify internal and external edges
  eext = (int*)malloc(numTriangleVertices*sizeof(int));
  partner = (int*)malloc(numTriangleVertices*sizeof(int));
  eext[edgelist[numTriangleVertices-1].t1] = 1; //Shortest edge is external by default
  partner[edgelist[numTriangleVertices-1].t1] = -1;
  for(vnum=0; vnum < numTriangleVertices-1; vnum++){
    if((edgelist[vnum].p1 == edgelist[vnum+1].p1) && (edgelist[vnum].p2 == edgelist[vnum+1].p2))
    {
      eext[edgelist[vnum].t1] = eext[edgelist[vnum+1].t1] = 0;
      partner[edgelist[vnum].t1] = edgelist[vnum+1].t1;
      partner[edgelist[vnum+1].t1] = edgelist[vnum].t1;
      vnum++;
    }
    else{
      eext[edgelist[vnum].t1] = 1;
      partner[edgelist[vnum].t1] = -1;
    }
  }
  
  //Remove external edges that can be removed
  removed = (int*)malloc(numTriangleVertices*sizeof(int));
  for(vnum=0;vnum<numTriangleVertices;vnum++) removed[vnum]=0;
  for(vnum=0;vnum<numTriangleVertices;vnum++){
    edge = edgelist[vnum].t1;
    edgeA = edge+1;
    if(edgeA%3 == 0) edgeA-=3;
    edgeB = edgeA+1;
    if(edgeB%3 == 0) edgeB-=3;
    if(!removed[edge] && len[edge]>minlen && eext[edge] &&
      TriangleCount[indexlist[edge]]>1 && TriangleCount[indexlist[edgeA]]>1 && TriangleCount[indexlist[edgeB]]>1){
      removed[edge] = removed[edgeA] = removed[edgeB] = 1;
      if(!eext[edgeA]) eext[partner[edgeA]] = 1;
      if(!eext[edgeB]) eext[partner[edgeB]] = 1;
      TriangleCount[indexlist[edge]]--;
      TriangleCount[indexlist[edgeA]]--;
      TriangleCount[indexlist[edgeB]]--;
      vnum = -1; //Very inefficient - returns to the beginning after each edge removal
    }
  }
  
  //Construct clockwise simple polygons
  numused = 0;
  boundary = (double*)malloc(2*numTriangleVertices*sizeof(double));
  for(*numpoly=0;*numpoly<MAXPOLY;(*numpoly)++){
      for(vnum=0; vnum<numTriangleVertices; vnum++){
          if(eext[vnum] && !removed[vnum]) break; //First external edge
      }
      if(vnum == numTriangleVertices) break;
      
      next = vnum;
      numboundarypoints[*numpoly] = 0;
      
      do{
          boundary[2*numused] = points[2*indexlist[next]];
          boundary[2*numused+1] = points[2*indexlist[next]+1];
          removed[next] = 1;
          (numboundarypoints[*numpoly])++;
          numused++;
          next++;
          if(next%3 == 0) next -= 3;
          
          while(!eext[next]){
              next = partner[next];
              next++;
              if(next%3 == 0) next -= 3;
          }
          
          if(numused == numTriangleVertices){
              printf("Too many external points in boundarypoints.\n");
              exit(1);
          }
      }
      while(next!=vnum);
  }

  free(indexlist);
  free(edgelist);
  free(len);
  free(eext);
  free(removed);
  free(partner);
  free(TriangleCount);

  return boundary;
}


int doublecompare(const void *e1,const void *e2)
{
  double *x1 = (double*)e1;
  double *x2 = (double*)e2;
  
  if(*x1 < *x2) return -1;
  else if(*x1 > *x2) return 1;
  return 0;
}


void updategrid(double *boundary,int numboundarypoints,Array2D<long> &grid,int xstart,int ystart,int gridsize)
{
  int bottomrow,toprow;
  
  //Get lists of crossovers
  Array2D<double> crossovers(grid.getRowNum(),numboundarypoints);
  int *numcrossovers = (int*)malloc(grid.getRowNum()*sizeof(int));
  for(int i=0;i<grid.getRowNum();i++) numcrossovers[i] = 0;

  for(int i=0; i<numboundarypoints;i++){
    double x1 = boundary[2*i];
    double y1 = boundary[2*i+1];
    double x2,y2;
    if(i == numboundarypoints-1){
      x2 = boundary[0];
      y2 = boundary[1];
    }
    else{
      x2 = boundary[2*(i+1)];
      y2 = boundary[2*(i+1)+1];
    }
    if(y1==y2) continue;

    if(y1<y2){
      bottomrow = (int)ceil((y1 - ystart)/gridsize);
      if(ystart + (bottomrow-1)*gridsize >= y1) bottomrow--;
      toprow = (int)floor((y2 - ystart)/gridsize);
      if(ystart + toprow*gridsize == y2) toprow--;
      else if(ystart + (toprow+1)*gridsize < y2) toprow++;
    }
    else{
      bottomrow = (int)ceil((y2 - ystart)/gridsize);
      if(ystart + (bottomrow-1)*gridsize >= y2) bottomrow--;
      toprow = (int)floor((y1 - ystart)/gridsize);
      if(ystart + toprow*gridsize == y1) toprow--;
      else if(ystart + (toprow+1)*gridsize < y1) toprow++;
    }

    double m = (x2 - x1) / (y2 - y1);
    for(int j=bottomrow;j<=toprow;j++){
      double y = ystart + j*gridsize;
      crossovers(j,numcrossovers[j]) = x1 + m*(y - y1);
      numcrossovers[j]++;
    }
  }
  
  //Get bounding box of boundary
  double left,right,top,bottom;
  left = right = boundary[0];
  top = bottom = boundary[1];
  for(int i=1;i<numboundarypoints;i++){
    if(boundary[2*i] < left) left = boundary[2*i];
    else if(boundary[2*i] > right) right = boundary[2*i];
    if(boundary[2*i+1] < bottom) bottom = boundary[2*i+1];
    else if(boundary[2*i+1] > top) top = boundary[2*i+1];
  }
  
  //Sort crossovers
  bottomrow = (int)floor((bottom - ystart)/gridsize);
  toprow = (int)ceil((top - ystart)/gridsize);
  for(int i=bottomrow;i<=toprow;i++){
    if(numcrossovers[i]>0) qsort(crossovers(i), numcrossovers[i], sizeof(double), &doublecompare);
  }
  
  //Increment grid
  int leftcol = (int)ceilf((left - xstart)/gridsize);
  int rightcol = (int)floorf((right - xstart)/gridsize);
  for(int i=bottomrow;i<=toprow;i++){
    int count = 0;
    for(int j=leftcol;j<=rightcol;j++){
      double x = xstart + j*gridsize;
      while(count<numcrossovers[i] && (x > crossovers(i,count) || (x == crossovers(i,count) && count%2 == 0))) count++;
      if(count%2 == 1) grid(j,i)++;
    }
  }
  
  free(numcrossovers);
}


int readsample(FILE *nestCoords,float *currentList)
{
    char line[1000];
    int numpts;
    
    
    fgets(line,1000,nestCoords);
    if(line[0] != 'P'){
        printf("Input file not at start of point set in readsample.\n");
        exit(1);
    }
    
    numpts = 0;
    while(fgets(line,1000,nestCoords) && strlen(line) > 1){
        sscanf(line,"%f %f",&currentList[2*numpts],&currentList[2*numpts+1]);
        currentList[2*numpts] -= XOFF;
        currentList[2*numpts+1] -= YOFF;
        numpts++;
    }
    
    return numpts;
}


void simpleBoundaries(FILE *nestCoords){
    double *boundary;
    int numboundarypoints[MAXPOLY],numpoly;
    int sample;//counter to iterate through the (non burn in) samples;
    int xstart = (MINX/GRIDSIZE)*GRIDSIZE;
    if(xstart < MINX) xstart=(MINX/GRIDSIZE+1)*GRIDSIZE;
    int ystart = (MINY/GRIDSIZE)*GRIDSIZE;
    if(ystart < MINY) ystart = (MINY/GRIDSIZE+1)*GRIDSIZE;
    int xwidth = (MAXX - xstart + 1)/GRIDSIZE;
    int ywidth = (MAXY - ystart + 1)/GRIDSIZE;
    xstart -= XOFF;
    ystart -= YOFF;
  
    Array2D<long> grid(xwidth,ywidth);
    for(int i=0;i<xwidth;i++){
      for(int j=0;j<ywidth;j++){
        grid(i,j) = 0;
      }
    }
  
    for (sample = 0; sample < SAMPLESIZE; sample++){
        float currentList[MAXNESTS];
        int numnests = readsample(nestCoords, currentList);

        printf("Processing sample %d. numnests=%d\n",sample+1,numnests);

        boundary = boundarypoints(currentList,numnests,numboundarypoints,&numpoly,MINLEN);
        
        printf("numpoly=%d\n",numpoly);
      
        int sum = 0;
        for(int i=0;i<numpoly;i++){
            updategrid(boundary+sum,numboundarypoints[i],grid,xstart,ystart,GRIDSIZE);
            sum += 2 * numboundarypoints[i];
        }
        free(boundary);
    }

    
    //Output heat map
    FILE *heat;
    heat = fopen("heatmap","w");
    if(!heat){
        printf("Could not open file heatmap in simpleBoundaries.\n");
        exit(1);
    }
    
    for(int i=0;i<xwidth;i++){
        for(int j=0;j<ywidth;j++){
            fprintf(heat,"%ld ",grid(i,j));
        }
        fprintf(heat,"\n");
    }
    
    fclose(heat);
    
    
    //Create polygons
    double threshvals[7] = {0.999,0.99,0.975,0.75,0.5,0.25,0.025};
    float *list = (float*)malloc(2*xwidth*ywidth*sizeof(float));
  
    //Output
    FILE* output = fopen(outputFileName,"w");
    if (!output){
        printf("output file failed to open.\n");
        exit(1);
    }

//    fprintf(output, "square = Polygon[{{%d,%d},{%d,%d},{%d,%d},{%d,%d}}];\n",
//            MINX, MINY, MAXX, MINY, MAXX, MAXY, MINX, MAXY);

    for(int t=0;t<7;t++){
        int numpoints = 0;
        int thresh = (int)floor((1.0-threshvals[t])*SAMPLESIZE);
        for(int i=0;i<xwidth;i++){
            for(int j=0;j<ywidth;j++){
                if((grid(i,j) > thresh) && (i==0 || i==xwidth-1 || j==0 || j==ywidth-1 ||
                                            (grid(i-1,j) <= thresh) || (grid(i+1,j) <= thresh) ||
                                            (grid(i,j-1) <= thresh) || (grid(i,j+1) <= thresh))){
                    list[2*numpoints] = (float)(xstart + i*GRIDSIZE);
                    list[2*numpoints+1] = (float)(ystart + j*GRIDSIZE);
                    numpoints++;
                }
            }
        }
        
        boundary = boundarypoints(list,numpoints,numboundarypoints,&numpoly,MINLEN);

        int k=0;
        for(int j=0;j<numpoly;j++){
            //fprintf(output, "%s-p%dcoord.%d = Polygon[{",nestCoordsFileName,(int)(1000*threshvals[t]),j+1);
            fprintf(output, "#%s-p%dcoord.%d\n x <- c(",nestCoordsFileName,(int)(1000*threshvals[t]),j+1);
            for (int i = 0; i < numboundarypoints[j]; i++){
                //fprintf(output, "{%f,%f}", boundary[2*k] + XOFF,boundary[2*k+1] + YOFF);
                fprintf(output, "%f,%f", boundary[2*k] + XOFF,boundary[2*k+1] + YOFF);
                if (i != numboundarypoints[j] - 1){
                    fprintf(output, ",");
                }
                if(i%50 == 49){
                    fprintf(output,"\n");
                }
                k++;
            }
            //fprintf(output, "}];\n");
            fprintf(output, ")\n");
        }

        free(boundary);
    }

    free(list);
    fclose(output);
}


int main(int argc, char **argv) {
    handleargs(argc, argv);

    FILE *nestCoords;//reading in the nest Coordinates
    nestCoords = fopen(nestCoordsFileName, "r");
    if (!nestCoords){
        printf("Could not open file %s.\n", nestCoordsFileName);
        exit(1);
    }
    else{
        printf("Opened file %s.\n", nestCoordsFileName);
    }

    simpleBoundaries(nestCoords);
    
    fclose(nestCoords);
}
