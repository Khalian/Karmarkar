//main problem with this program as of now is the input data type(int or float) and associated type cast errors.
//we ultimately decided to keep the function void, we have absolutely no idea yet about the exact parameters or arguments to pass
//need to make d, get input of constraints and most importantly change the way this program displays output
//a=m*n matrix

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
 
uint32_t         stampstart();
uint32_t         stampstop(uint32_t start);

void phi(float c[50][50],float a[50][50],float d[50][50],float x[50][1],int n,int m);
float minorm(float arr[50][50],int row,int col,int order);                  
void transpose(float a[50][50],int n,int m);
void inverse(float a[50][50],int n);
float determinant(float b[50][50],int m);
void subtraction(float a[50][50],float b[50][50],int row,int col);
float mod(float arr[50][50],int row,int col);
void mul(float a[50][50],float b[50][50],int n,int m, int p,int q);
void augment(float arr[50][50],int m,int n);
int compvals(float x1[50][1],float x2[50][1],int n);

float alpha,delta;

float f(float x[50][1],float c[50][50],int n); 
float consradii(int n);
void display(float c[50][50],float x[50][1],int n); 

int main()                        //this function actually runs the karmarkar's algorithm, ie minimising objective funtion
{
    int n,m,i,j,k,L,ceil1,ceil2,K,q,cmp;
    float x1[50][1],x2[50][1],x0[50][1],P,c[50][50],c2[50][50],temp[50][50],log2,logP,logn,val,f1,f2;
    float valx,val0;
    float a[50][50],d[50][50],aintact[50][50];
    uint32_t         start, stop;

    printf("\n Enter the order of the linear program:");
    scanf("%d",&n);

    printf("\n Enter the number of constraints:");
    scanf("%d",&m);

    printf("\n Enter the coefficients of the variables in the objective funtion:");
    for(i=0;i<n;i++)                  //input of objective function
        scanf("%f",&c[i][0]);             //c actually behaves like a row matrix, inspite of its definition         

    printf("\n Enter the matrix of constraints rowwise:");
    for(i=0;i<m;i++)                 //loop for input of constarints
        for(j=0;j<n;j++)
        {
            scanf("%f",&a[i][j]); 
            aintact[i][j]=a[i][j];
        }

    start = stampstart();

    alpha=n-1;                        //initialising alpha
    alpha=(alpha)/(3*n);

    delta=alpha*alpha;                //calculation of delta
    delta=delta/(1-alpha);
    delta=alpha-delta;

    L=n*m;

    P=1.0;

    for(i=0;i<n;i++)                  //calculation of P,step needs review
    {
        if(c[i][0]!=0)
            P=P*c[i][0];
    }

    P=abs(P);

    log2=log(2);
    logP=log(P)/log2;
    logn=log(n)/log2;
    ceil1=(int)(logP);
    ceil1=ceil1+1;

    ceil2=(int)(logn);
    ceil2=ceil2+1;
    ceil2=ceil2*n;

    L=L+ceil1+ceil2;                  //L is an int, calculation of L
    K=(int)((2*n*L)/(delta));
    K++;                              // iteration bound K calculated

    val=(float)(n);
    val=1/val;                          //i dont know if it has type cast errors

    for(i=0;i<n;i++)                  //initialisation of x0
        x0[i][0]=val;
    
    for(i=0;i<n;i++)                 // c shall always be intact throughout the program,c2 is the temp c
        c2[i][0]=c[i][0];

    for(i=0;i<n;i++)                 //for compatibility in mul fn
        temp[i][0]=x0[i][0];

    transpose(c2,n,1);               //transposing the temporary c

    mul(c2,temp,1,n,n,1);              // the base case

    val=c2[0][0];

    if(val==0)                         //check if a0=x0 is the answer
    {
        display(c,x0,n);
        stop=stampstop(start);
        exit(0);
    }

    for(i=0;i<n;i++)                  //this loop initialises the first x1 as x0, ie the starting point
        x1[i][0]=x0[i][0];

    for(i=0;i<K;i++)                 //this is the actual loop to calculate the optimal point
    {

        for(j=0;j<n;j++)                 //x2 will now posses values of x1, x2 goes under transformation
        x2[j][0]=x1[j][0];

        for(j=0;j<n;j++)
        {
            for(k=0;k<n;k++)
            {
                if(k==j)
                    d[j][k]=x2[j][0];
                else
                    d[j][k]=0;
            }
        }

        phi(c,aintact,d,x2,n,m);               //here i am phi ing x2 while keeping the original point as x1 

        for(j=0;j<n;j++)                 // c shall always be intact throughout the program
             c2[j][0]=c[j][0];

        transpose(c2,n,1);               //transposing the temporary c

        for(j=0;j<n;j++)                 //because x2 is a column matrix and hence is unsuitable for multiplication
            temp[j][0]=x2[j][0];

        mul(c2,temp,1,n,n,1);

        if(val==0)
        {         
            display(c,x2,n);
            stop=stampstop(start);
            exit(0);
        }

        cmp=compvals(x1,x2,n);
        if(cmp==0)
        {
            display(c,x2,n);
            stop=stampstop(start);
            exit(0);
        }

        f2=f(x2,c,n);                    //the positive optimality test
        f1=f(x1,c,n);
        /*
        printf("\n x1 is:");
        for(j=0;j<n;j++)
        printf(" %f",x1[j][0]);

        printf("\n x2 is:");
        for(j=0;j<n;j++)
        printf(" %f",x2[j][0]);

        printf("\n Values of fx, fx-1 and delta  are resp %f %f %f",f2,f1,delta);
        exit(0);
        */

        if(f2>(f1-(0.01*delta)))                //needs reviewing
        {
            printf("\n The solution is unfeasible or unbounded, no finite optimal solution exists for it");
            stop=stampstop(start);
            exit(0);
        }
                                  //under repair

    //till here we have done as required till step2 of the linear programming book

        for(j=0;j<n;j++)
            x1[j][0]=x2[j][0];             //for continuance of the loop ie x1 becomes the generated x2, step needs review

        //the termination condition q as defined in the paper, i dont know what value it takes but i think its 2*L

        for(j=0;j<n;j++)                 // c shall always be intact throughout the program,c2 is the temp c
            c2[j][0]=c[j][0];

        transpose(c2,n,1);

        for(j=0;j<n;j++)
            temp[j][0]=x2[j][0];

        mul(c2,temp,1,n,n,1);
        valx=c2[0][0];

        for(j=0;j<n;j++)                 // c shall always be intact throughout the program,c2 is the temp c
            c2[j][0]=c[j][0];

        transpose(c2,n,1);

        for(j=0;j<n;j++)
            temp[j][0]=x0[j][0];

        mul(c2,temp,1,n,n,1);
        val0=c2[0][0];

        valx=valx/val0;

    }                               // end of principal loop

    display(c,x2,n);
    stop=stampstop(start);
}                                //end of main
 
void phi(float c[50][50],float a[50][50],float d[50][50],float x[50][1],int n,int m) // points are stored by x
                                                                          //this function calculates the transformation of
{                                                                         //karmarkar's algorithm
    float cdash[50][50];                //need serious reviewing                                                    
    float cp[50][50];
    float z[50][50];
    float cpcap[50][50];
    int i,j;
    float value;
    float anot[50][50];
    float bprime[50][50];
    float e[50][50];
    float ddash[50][50];
    float constant;
    float b[50][50];
    float val,atemp[50][50];

    for(i=0;i<m;i++)                 //loop for input of constarints
        for(j=0;j<n;j++)
            atemp[i][j]=a[i][j];

    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            ddash[i][j]=d[i][j];             //copying d to ddash since d will b destroyed

    mul(d,c,n,n,n,1);                //d=d.c

    for(i=0;i<n;i++)
        cdash[i][0]=d[i][0];             //c'=d's first column

    mul(atemp,ddash,m,n,n,n);            //a=a.d 

    augment(atemp,m,n);                  //a=augment of a

    for(i=0;i<(m+1);i++)                //z is a m+1*n matrix so is a   
        for(j=0;j<n;j++)
            z[i][j]=atemp[i][j];                  //z=b

    transpose(atemp,m+1,n);               //atemp=bt

    mul(z,atemp,m+1,n,n,m+1);         //z=zztrans=BBtrans 


    inverse(z,m+1);               //z=zinverse

    transpose(atemp,n,m+1);           //taking the transpose of a ie Btrans to obtain back the original z or B ie a=B

    mul(z,atemp,m+1,m+1,m+1,n);       //z=inv(zztrans)z or z=inv(BBtrans)B

    mul(z,cdash,m+1,n,n,1);       //multiplying with cdash result in z

    transpose(atemp,m+1,n);           //taking transpose of a ie z or Btrans

    mul(atemp,z,n,m+1,m+1,1);         //m loses its effects from here onwards  


    subtraction(cdash,atemp,n,1);     //subtracting a from cdash result in cdash
    
    for(i=0;i<n;i++)
        cp[i][0]=cdash[i][0];         //copying cdash to  cp cp=n*1 matrix

    value=mod(cp,n,1);            //obtaining the mod of cp

    if(value==0)
    {
        for(i=0;i<n;i++)
        cpcap[i][0]=0;         //cpcap=0;for all zeros
    }

    else
    {
        value=1/value;
        for(i=0;i<n;i++)
        cpcap[i][0]=cp[i][0]*value;
    }                              //obtaing cpcap from formula, cpcap=n*1 matrix

    value=(float)n;
    value=1/value;

    for(i=0;i<n;i++)
        anot[i][0]=value;               //initializing anot to 1\n

    val=consradii(n);
    val=alpha*val;

    for(i=0;i<n;i++)
        cpcap[i][0]=(val*cpcap[i][0]);

    subtraction(anot,cpcap,n,1);  //subtracting anot from cpcap result in anot

    for(i=0;i<n;i++)
       bprime[i][0]=anot[i][0];      //copying anot to bprime

    for(i=0;i<n;i++)
        e[0][i]=1;                    //initializing the elements of the array e with 1

    mul(e,ddash,1,n,n,n);         //multiplying e with ddash or the original d
    mul(e,bprime,1,n,n,1);
    constant=e[0][0];             //the single element array in e is terrmed constant

    mul(ddash,bprime,n,n,n,1);

    for(i=0;i<n;i++)
    {
        b[i][0]=ddash[i][0];          //assigns the values of ddash to b
        b[i][0]=b[i][0]/constant;     //computes the value of b
    }

    for(i=0;i<n;i++)
        x[i][0]=b[i][0];
}

void inverse(float arr[50][50],int n)        // job done, cross checked
{
    int i,j;
    float fdet;
    float inv[50][50];

    if(n==2)
    {
        fdet=(arr[0][0]*arr[1][1])-(arr[0][1]*arr[1][0]);
        fdet=1/fdet;
        inv[0][0]=(arr[1][1])*fdet;
        inv[0][1]=(-arr[1][0])*fdet;
        inv[1][0]=(-arr[0][1])*fdet;
        inv[1][1]=(arr[0][0])*fdet;

    for(i=0;i<2;i++)                            //final assignment loop
        for(j=0;j<2;j++)
            arr[i][j]=inv[i][j];
    goto end;
}

    fdet=determinant(arr,n);

    if(fdet==0)
    {printf("\n Matrix is singular, no inverse exists");
    exit(0);}

    fdet=1/fdet;

    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            inv[i][j]=minorm(arr,i,j,n);
            if((i+j)%2!=0)
            inv[i][j]=-inv[i][j];
            inv[i][j]=inv[i][j]*(fdet);
        }

    for(i=0;i<n;i++)                            //final assignment loop
        for(j=0;j<n;j++)
            arr[i][j]=inv[i][j];
            end:
            printf("");
}

float minorm(float arr[50][50],int row,int col,int order)     //job done,cross checked
{
    float result;
    int i,j,a,b;
    float auxarr[50][50];
    order=order-1;         
    i=j=0;

    for(a=0;a<order;a++)
    {
        j=0;
        if(i==row)
            i++;
        for(b=0;b<order;b++)
        {
            if(j==col)
            j++;
            auxarr[a][b]=arr[i][j];
            j++;
        }               //end of b loop 
            i++;
        }               //end of a loop
        result=determinant(auxarr,order);
    return(result);
}

float determinant(float b[50][50],int m)                       //job done, cross checked 
{
    int i,j,p,sign=1,k;
    float c[50][50];
    float sum=0;
    if(m==2)
    {
        sum = b[0][0]*b[1][1] - b[0][1]*b[1][0];
        return sum;
    }
    for(p=0;p<m;p++)
    {
        int h = 0,k = 0;
        for(i=1;i<m;i++)
    {
    for( j=0;j<m;j++)
    {
        if(j==p)
        continue;
        c[h][k] = b[i][j];
        k++;
        if(k == m-1)
        {
            h++;
            k = 0;
        }
    }
}
    for(k=0;k<p;k++)
    sign=sign*(-1);

    sum=sum+b[0][p]*sign*determinant(c,m-1);
}
    return sum;
}

void transpose(float a[50][50],int n,int m)                 //a[][] is transformed care when using it in your phi function
{
    int i,j;
    float b[50][50];
    for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    b[i][j]=a[j][i];

    for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    a[i][j]=0;

    for(i=0;i<m;i++)
    for(j=0;j<n;j++)
    a[i][j]=b[i][j];
}

void subtraction(float a[50][50],float b[50][50],int row,int col)            //converts a[][]=a[][]-b[][] care when using
{
    int i,j;
    for(i=0;i<row;i++)
    for(j=0;j<col;j++)
    a[i][j]=a[i][j]-b[i][j];
}

void augment(float arr[50][50],int m,int n)                          //m no of rows n no of columns
{
    int i;
    for(i=0;i<n;i++)
    arr[m][i]=1;
}

float mod(float arr[50][50],int row,int col)
{
    float arr2[50][50];
    float res;
    res=0;
    int i,j;
    for(i=0;i<row;i++)
        for(j=0;j<col;j++)
        {
            arr2[i][j]=arr[i][j]*arr[i][j];
            res=res+arr2[i][j];
        }

    res=sqrt(res);
    return(res);
}

void mul(float a[50][50],float b[50][50],int n,int m, int p,int q)    //mul is replacing a with a*b
{
    float c[50][50];
    int i,j,k;

    if(p!=m)
    {
        printf("\n Matrix multiplication not possible");
        exit(0);
    }

    for(i=0;i<n;i++)    
        for(j=0;j<q;j++)
            c[i][j]=0;

    for(i=0;i<n;i++)
        for(j=0;j<q;j++)
            for(k=0;k<m;k++)
                c[i][j]=c[i][j]+a[i][k]*b[k][j];

    for(i=0;i<n;i++)
        for(j=0;j<q;j++)
            a[i][j]=c[i][j];
        }

float f(float x[50][1],float c[50][50],int n)           //need to cross check
{
    float result,temp[50][50],ctemp[50][50];
    float var=0,prod=1;
    int i;

    for(i=0;i<n;i++)
    ctemp[i][0]=c[i][0];

    transpose(ctemp,n,1);                      //c becomes c transpose

    for(i=0;i<n;i++)
    temp[i][0]=x[i][0];                    //done for making multiplication possible,check the prototype of the mul fn

    mul(ctemp,temp,1,n,n,1);                      //c=ct*x
    var=ctemp[0][0];                           //var=element of ctx which is a 1*1 matrix

    for(i=0;i<n;i++)
    {
        if(x[i][0]!=0)                         //doubtful step
            prod=prod*x[i][0];
    }
    result=var/prod;
    result=log(result);
    result=n*result;

    return(result);

} 

float consradii(int n)                  //tested
{
    float r;
    r=(n*(n-1));
    r=1/r;
    r=sqrt(r);
    return r;
}

void display(float c[50][50],float x[50][1],int n)              //not tested
{
    float temp[50][50],result,c2[50][50];
    int i;

    for(i=0;i<n;i++)                                          //for the mul func
        temp[i][0]=x[i][0];

    for(i=0;i<n;i++)
        c2[i][0]=c[i][0];

    transpose(c2,n,1);
    mul(c2,temp,1,n,n,1);

    result=c2[0][0];

    printf("\n The point at which optimality is obtained is:");

    for(i=0;i<n;i++)
    printf(" %f",x[i][0]);

    printf("\n Well folks it has been a wonderful programming marathon, neways, z=cTx has the minimum value %f",result);
}

int compvals(float x1[50][1],float x2[50][1],int n)
{
int flag=0;
int i;
for(i=0;i<n;i++)
{
if((x1[i][0]-x2[i][0])>0.001)
flag=1;
}
return(flag);
}

uint32_t
stampstart()
{
	struct timeval  tv;
	struct timezone tz;
	struct tm      *tm;
	uint32_t         start;
 
	gettimeofday(&tv, &tz);
	tm = localtime(&tv.tv_sec);
 
	printf("TIMESTAMP-START\t  %d:%02d:%02d:%d (~%d ms)\n", tm->tm_hour,
	       tm->tm_min, tm->tm_sec, tv.tv_usec,
	       tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
	       tm->tm_sec * 1000 + tv.tv_usec / 1000);
 
	start = tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
		tm->tm_sec * 1000 + tv.tv_usec / 1000;
 
	return (start);
 
}
 
uint32_t
stampstop(uint32_t start)
{
 
	struct timeval  tv;
	struct timezone tz;
	struct tm      *tm;
	uint32_t         stop;
 
	gettimeofday(&tv, &tz);
	tm = localtime(&tv.tv_sec);
 
	stop = tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
		tm->tm_sec * 1000 + tv.tv_usec / 1000;
 
	printf("TIMESTAMP-END\t  %d:%02d:%02d:%d (~%d ms) \n", tm->tm_hour,
	       tm->tm_min, tm->tm_sec, tv.tv_usec,
	       tm->tm_hour * 3600 * 1000 + tm->tm_min * 60 * 1000 +
	       tm->tm_sec * 1000 + tv.tv_usec / 1000);
 
	printf("ELAPSED\t  %d ms\n", stop - start);
 
	return (stop);
 
}









